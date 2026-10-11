#!/usr/bin/env python
"""
Pure-Python tabix client over HTTP range requests.
No pysam and no tabix binary (neither is installable on this Windows host).

Implements from scratch:
  - BGZF block detection + decompression (gzip member carrying the 'BC' extra subfield)
  - .tbi index parsing (htslib tabix index spec)
  - reg2bins binning scheme + linear index lookup
  - region query served purely by HTTP Range, never downloading the whole data file

Every query returns (rows, stats); stats records exactly how many bytes crossed
the wire and what fraction of the remote file that was, so cost claims are auditable.
"""
import gzip
import json
import struct
import sys
import time
import zlib

import requests

MAX_BGZF_BLOCK = 65536


# ---------------------------- BGZF ----------------------------
def is_bgzf(first_bytes):
    """Return (bool, bsize) given at least the first 18 bytes of a gzip member."""
    b = first_bytes
    if len(b) < 18 or b[0] != 0x1F or b[1] != 0x8B or b[2] != 8:
        return False, None
    flg = b[3]
    if not (flg & 0x04):  # FEXTRA is mandatory for BGZF
        return False, None
    xlen = struct.unpack("<H", b[10:12])[0]
    extra = b[12:12 + xlen]
    i = 0
    while i + 4 <= len(extra):
        si1 = extra[i]
        si2 = extra[i + 1]
        slen = struct.unpack("<H", extra[i + 2:i + 4])[0]
        if si1 == 66 and si2 == 67 and slen == 2:  # ord('B'), ord('C')
            bsize = struct.unpack("<H", extra[i + 4:i + 6])[0] + 1
            return True, bsize
        i += 4 + slen
    return False, None


def _decompress_consumed(buf):
    """Decompress whole BGZF blocks; return (data, bytes_consumed).

    bytes_consumed counts only complete blocks, so a sliced read resumes exactly at a
    block boundary.
    """
    out = bytearray()
    pos = 0
    while pos + 18 <= len(buf):
        ok, bsize = is_bgzf(buf[pos:pos + 18])
        if not ok or pos + bsize > len(buf):
            break
        out += zlib.decompress(buf[pos:pos + bsize], 16 + zlib.MAX_WBITS)
        pos += bsize
    return bytes(out), pos


def bgzf_decompress_stream(buf):
    """Decompress consecutive BGZF blocks from buf. Returns (data, n_blocks)."""
    out = bytearray()
    pos = 0
    nblk = 0
    while pos + 18 <= len(buf):
        ok, bsize = is_bgzf(buf[pos:pos + 18])
        if not ok or pos + bsize > len(buf):
            break
        out += zlib.decompress(buf[pos:pos + bsize], 16 + zlib.MAX_WBITS)
        pos += bsize
        nblk += 1
    return bytes(out), nblk


# ------------------------- tabix index -------------------------
class Tbi(object):
    def __init__(self, raw):
        d = gzip.decompress(raw)
        if d[:4] != b"TBI\x01":
            raise ValueError("not a TBI index, magic=%r" % d[:4])
        (self.n_ref, self.fmt, self.col_seq, self.col_beg, self.col_end,
         self.meta, self.skip, l_nm) = struct.unpack("<8i", d[4:36])
        p = 36
        self.names = [n.decode() for n in d[p:p + l_nm].split(b"\x00") if n]
        p += l_nm
        self.bins = []
        self.linear = []
        for _ in range(self.n_ref):
            n_bin = struct.unpack("<i", d[p:p + 4])[0]
            p += 4
            bd = {}
            for _ in range(n_bin):
                b, n_chunk = struct.unpack("<Ii", d[p:p + 8])
                p += 8
                ch = []
                for _ in range(n_chunk):
                    beg, end = struct.unpack("<QQ", d[p:p + 16])
                    p += 16
                    ch.append((beg, end))
                bd[b] = ch
            n_intv = struct.unpack("<i", d[p:p + 4])[0]
            p += 4
            lin = list(struct.unpack("<%dQ" % n_intv, d[p:p + 8 * n_intv])) if n_intv else []
            p += 8 * n_intv
            self.bins.append(bd)
            self.linear.append(lin)

    @staticmethod
    def reg2bins(beg, end):
        """0-based half-open interval -> candidate bin ids (htslib scheme)."""
        lst = [0]
        end -= 1
        for shift, off in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
            lst += list(range(off + (beg >> shift), off + (end >> shift) + 1))
        return lst

    def chunks(self, chrom, beg, end):
        """beg/end 0-based half-open. Returns merged (start, end) virtual offsets.

        Both a lower and an UPPER bound come from the linear index. The upper
        bound matters: coarse bins (especially bin 0, which VCFs populate with
        chromosome-spanning SV records) carry chunks that run to the end of the
        contig, so taking max(chunk_end) alone would ask the server for hundreds
        of MB to answer a 200 kb query.
        """
        if chrom not in self.names:
            return []
        i = self.names.index(chrom)
        lin = self.linear[i]
        li = beg >> 14
        if li < len(lin):
            min_off = lin[li]
        elif lin:
            min_off = lin[-1]
        else:
            min_off = 0
        # The linear index is a valid LOWER bound only. Using lin[] as an upper bound
        # is unsound: a long record (e.g. a structural variant) that overlaps a later
        # 16 kb interval drags that interval's entry back to an early file offset, so
        # the "upper bound" can fall below data that genuinely belongs to the query.
        # That truncated regions silently -- 0 of 1226 eQTLGen hit positions were
        # recovered for SLC2A9 before this was removed. Callers must fetch each chunk
        # and filter by position, which is what tabix itself does.
        cand = []
        for b in self.reg2bins(beg, end):
            for cb, ce in self.bins[i].get(b, []):
                lo = max(cb, min_off)
                if ce > lo:
                    cand.append((lo, ce))
        cand.sort()
        out = []
        for c in cand:
            if out and c[0] <= out[-1][1]:
                out[-1][1] = max(out[-1][1], c[1])
            else:
                out.append([c[0], c[1]])
        return [tuple(x) for x in out]


# ---------------------------- query ----------------------------
def http_get_range(url, start, end, sess=None, stats=None, retries=4):
    """Ranged GET, read as a stream. The 1000genomes host drops long-running
    connections mid-body, so short reads are retried rather than silently
    returning a truncated buffer."""
    s = sess or requests
    want = end - start + 1
    last = None
    for attempt in range(retries):
        try:
            r = s.get(url, headers={"Range": "bytes=%d-%d" % (start, end)},
                      timeout=300, stream=True)
            if r.status_code not in (200, 206):
                raise RuntimeError("HTTP %d for %s [%d-%d]"
                                   % (r.status_code, url, start, end))
            buf = bytearray()
            for chunk in r.iter_content(1 << 20):
                buf += chunk
            body = bytes(buf)
            if stats is not None:
                stats["requests"] = stats.get("requests", 0) + 1
                stats["bytes"] = stats.get("bytes", 0) + len(body)
                stats.setdefault("status_codes", []).append(r.status_code)
                if attempt:
                    stats["retries"] = stats.get("retries", 0) + attempt
            # a short read on a 206 means the connection died mid-body
            if r.status_code == 206 and len(body) < want and attempt < retries - 1:
                last = RuntimeError("short read %d/%d" % (len(body), want))
                time.sleep(2 * (attempt + 1))
                continue
            return body, r
        except (requests.exceptions.ChunkedEncodingError,
                requests.exceptions.ConnectionError,
                requests.exceptions.ReadTimeout) as e:
            last = e
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("range GET failed after %d attempts: %r" % (retries, last))


def query(url, chrom, start1, end1, tbi_url=None, sess=None):
    """1-based inclusive region query. Returns (rows, stats)."""
    sess = sess or requests.Session()
    stats = {"bytes": 0, "requests": 0, "data_url": url}
    tbi_url = tbi_url or url + ".tbi"
    stats["tbi_url"] = tbi_url

    # 0) does the server honour byte ranges on the DATA file?
    head = sess.head(url, timeout=60)
    stats["head_status"] = head.status_code
    stats["accept_ranges"] = head.headers.get("Accept-Ranges")
    stats["data_file_total_bytes"] = int(head.headers.get("Content-Length", -1))

    # 1) is the data file truly BGZF (not plain gzip)? read 18 bytes only.
    probe, r = http_get_range(url, 0, 17, sess, stats)
    ok, bsize = is_bgzf(probe)
    stats["data_is_bgzf"] = ok
    stats["first_block_size"] = bsize
    stats["probe_status"] = r.status_code
    if not ok:
        raise RuntimeError("data file is NOT BGZF; tabix range query is impossible")

    # 2) pull the index (small) and parse it
    ri = sess.get(tbi_url, timeout=300)
    if ri.status_code != 200:
        raise RuntimeError("HTTP %d fetching index %s" % (ri.status_code, tbi_url))
    stats["requests"] += 1
    stats["bytes"] += len(ri.content)
    stats["tbi_bytes"] = len(ri.content)
    tbi = Tbi(ri.content)
    stats["tbi_n_ref"] = tbi.n_ref
    stats["tbi_names_head"] = tbi.names[:5]
    stats["tbi_cols_seq_beg_end"] = [tbi.col_seq, tbi.col_beg, tbi.col_end]

    # 3) index lookup -> byte span
    use = chrom
    if use not in tbi.names:
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alt in tbi.names:
            use = alt
        else:
            raise RuntimeError("chrom %r absent from index; index names start %r"
                               % (chrom, tbi.names[:5]))
    stats["chrom_used"] = use
    ch = tbi.chunks(use, start1 - 1, end1)
    stats["n_chunks"] = len(ch)
    if not ch:
        stats["rows"] = 0
        return [], stats
    # 4) Fetch each chunk as its own sliced byte range, stopping once records pass the
    # region end. Collapsing disjoint chunks into one min->max span requested 491 MB for
    # a 230 kb window on VCFs whose windows are scattered by structural-variant records.
    GAP = 1 << 16
    spans = []
    for cb_, ce_ in sorted(ch):
        lo, hi = cb_ >> 16, (ce_ >> 16) + MAX_BGZF_BLOCK
        if spans and lo - spans[-1][1] <= GAP:
            spans[-1][1] = max(spans[-1][1], hi)
        else:
            spans.append([lo, hi])
    stats["byte_span"] = [spans[0][0], spans[-1][1]]
    stats["n_spans"] = len(spans)

    SLICE = 1 << 20
    parts = []
    nblk = 0
    for lo, hi in spans:
        cur = lo
        carry = b""
        while cur < hi:
            buf, r = http_get_range(url, cur, min(cur + SLICE - 1, hi), sess, stats)
            if not buf:
                break
            cdata, consumed = _decompress_consumed(buf)
            if consumed == 0:
                break
            # a record can straddle a BGZF block boundary; carry the partial line over
            cdata = carry + cdata
            lines = cdata.split(b"\n")
            carry = lines.pop()
            parts.extend(lines)
            cur += consumed
            nblk += 1
            past = False
            for l in reversed(lines):
                if not l or l[:1] == b"#":
                    continue
                f2 = l.split(b"\t", 3)
                if len(f2) < 3:
                    continue
                try:
                    past = int(f2[1]) > end1
                except ValueError:
                    continue
                break
            if past:
                break
    data = b"\n".join(parts)
    stats["blocks_decompressed"] = nblk
    stats["uncompressed_bytes"] = len(data)

    # 5) parse + filter to the requested region
    cs = tbi.col_seq - 1
    cb = tbi.col_beg - 1
    meta = chr(tbi.meta) if 0 < tbi.meta < 128 else "#"
    zero_based = bool(tbi.fmt & 0x10000)
    rows = []
    for line in data.split(b"\n"):
        if not line:
            continue
        txt = line.decode("utf-8", "replace").rstrip("\r")
        if txt[:1] == meta:
            continue
        f = txt.split("\t")
        if len(f) <= max(cs, cb) or f[cs] != use:
            continue
        try:
            pos = int(f[cb])
        except ValueError:
            continue
        p1 = pos + 1 if zero_based else pos
        if start1 <= p1 <= end1:
            rows.append(f)
    stats["rows"] = len(rows)
    if stats["data_file_total_bytes"] > 0:
        stats["pct_of_file_transferred"] = round(
            100.0 * stats["bytes"] / stats["data_file_total_bytes"], 6)
    return rows, stats


def header_lines(url, sess=None, nbytes=65536):
    """Read just the first BGZF block(s) and return the leading text lines."""
    sess = sess or requests.Session()
    stats = {"bytes": 0, "requests": 0}
    buf, _ = http_get_range(url, 0, nbytes - 1, sess, stats)
    data, _ = bgzf_decompress_stream(buf)
    return data.decode("utf-8", "replace").split("\n"), stats


if __name__ == "__main__":
    u, c, s, e = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    t = sys.argv[5] if len(sys.argv) > 5 else None
    rows, st = query(u, c, s, e, tbi_url=t)
    print(json.dumps(st, indent=2))
    print("--- first 5 rows ---")
    for row in rows[:5]:
        print("\t".join(row))
    print("--- total rows in region: %d ---" % len(rows))
