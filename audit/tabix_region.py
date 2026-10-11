"""AUDIT: minimal remote tabix reader (HTTP range requests), written for this audit.

Does not import anything from the project's src/. Reads the .tbi, finds the BGZF chunks that can
overlap chrom:start-end, range-requests them and returns the matching text lines.

CLI:  python audit/tabix_region.py <url.tsv.gz> <chrom> <start> <end>
"""
import gzip, struct, sys, urllib.request, zlib


def _get(url, lo=None, hi=None):
    req = urllib.request.Request(url, headers={"User-Agent": "c2d-audit/1.0"})
    if lo is not None:
        req.add_header("Range", f"bytes={lo}-{hi}")
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except Exception as e:  # retry transient network errors
            err = e
    raise err


def _reg2bins(beg, end):
    end -= 1
    bins = [0]
    for shift, off in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        bins.extend(range(off + (beg >> shift), off + (end >> shift) + 1))
    return bins


def read_index(url):
    raw = gzip.decompress(_get(url + ".tbi"))
    assert raw[:4] == b"TBI\1"
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", raw[4:36])
    names = raw[36:36 + l_nm].rstrip(b"\0").split(b"\0")
    p = 36 + l_nm
    idx = {}
    for name in names:
        (n_bin,) = struct.unpack("<i", raw[p:p + 4]); p += 4
        bins = {}
        for _ in range(n_bin):
            b, n_chunk = struct.unpack("<Ii", raw[p:p + 8]); p += 8
            bins[b] = [struct.unpack("<QQ", raw[p + 16 * i:p + 16 * i + 16]) for i in range(n_chunk)]
            p += 16 * n_chunk
        (n_intv,) = struct.unpack("<i", raw[p:p + 4]); p += 4
        lin = struct.unpack(f"<{n_intv}Q", raw[p:p + 8 * n_intv]); p += 8 * n_intv
        idx[name.decode()] = (bins, lin)
    return dict(idx=idx, col_seq=col_seq, col_beg=col_beg, col_end=col_end, fmt=fmt, meta=chr(meta), skip=skip)


def fetch(url, chrom, start, end, index=None):
    """Lines whose position (1-based, column col_beg) lies in [start, end]."""
    ix = index or read_index(url)
    if chrom not in ix["idx"]:
        return [], ix
    bins, lin = ix["idx"][chrom]
    min_off = lin[min(start >> 14, len(lin) - 1)] if lin else 0
    chunks = [c for b in _reg2bins(start - 1, end) for c in bins.get(b, []) if c[1] > min_off]
    if not chunks:
        return [], ix
    lo = min(c[0] for c in chunks); hi = max(c[1] for c in chunks)
    c_lo, u_lo = lo >> 16, lo & 0xFFFF
    c_hi = (hi >> 16) + 65536 + 100          # the last block is at most 64 KiB compressed
    blob = _get(url, c_lo, c_hi)
    out, p, first = [], 0, True
    buf = b""
    while p + 18 <= len(blob):
        if blob[p:p + 4] != b"\x1f\x8b\x08\x04":
            break
        xlen = struct.unpack("<H", blob[p + 10:p + 12])[0]
        bsize = struct.unpack("<H", blob[p + 16:p + 18])[0] + 1
        if p + bsize > len(blob):
            break
        data = zlib.decompress(blob[p + 12 + xlen:p + bsize - 8], -15)
        buf += data[u_lo:] if first else data
        first = False
        p += bsize
        if c_lo + p > (hi >> 16) + 1 and len(data) == 0:
            break
    cs, cb = ix["col_seq"] - 1, ix["col_beg"] - 1
    for line in buf.split(b"\n"):
        f = line.split(b"\t")
        if len(f) <= max(cs, cb) or f[cs].decode() != chrom:
            continue
        try:
            pos = int(f[cb])
        except ValueError:
            continue
        if start <= pos <= end:
            out.append(line.decode())
    return out, ix


def header(url):
    """Column names: first line of the first BGZF block."""
    blob = _get(url, 0, 70000)
    xlen = struct.unpack("<H", blob[10:12])[0]
    bsize = struct.unpack("<H", blob[16:18])[0] + 1
    return zlib.decompress(blob[12 + xlen:bsize - 8], -15).split(b"\n")[0].decode().split("\t")


if __name__ == "__main__":
    u, c, s, e = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    print("\t".join(header(u)))
    for l in fetch(u, c, s, e)[0]:
        print(l)
