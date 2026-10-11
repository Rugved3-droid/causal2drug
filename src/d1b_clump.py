#!/usr/bin/env python
"""Deliverable 1b step 4: LD-clump eQTLGen cis-eQTL hits per Minikel target.

Pre-registered rules applied unchanged (R1/R5, DECISION 6):
  window   = gene body +- 100 kb  (GRCh37, Ensembl release 87 GTF)
  p        < 5e-8                 (FROZEN by DECISION 6)
  clumping = greedy, lowest p first, drop r2 >= 0.01
  LD panel = 1000 Genomes phase 3 (GRCh37), EUR, 503 samples

Performance note: the .tbi index and the VCF #CHROM sample header are fetched ONCE per
chromosome and reused across every target on it. Re-fetching the 176 KB index per target
dominated runtime otherwise.
"""
import collections
import concurrent.futures as cf
import csv
import json
import os
import sys
import threading
import zlib
import time

import numpy as np
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabix_http as T  # noqa: E402

D = r"D:\causal2drug\deliverables\d1b"
PANEL = r"D:\causal2drug\data\raw\integrated_call_samples_v3.20130502.ALL.panel"
KG = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
      "ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
FLANK = 100_000
R2_CLUMP = 0.01
WORKERS = 6
CKPT = D + r"\clump_ckpt"

BYTES = collections.Counter()


def eur_samples():
    with open(PANEL) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return {r["sample"] for r in rd if r["super_pop"] == "EUR"}


def load_inputs():
    coords = {r["symbol"]: r for r in csv.DictReader(
        open(D + r"\gene_coords_grch37.tsv", encoding="utf-8"), delimiter="\t")}
    hits = collections.defaultdict(list)
    with open(D + r"\eqtlgen_hits_targets.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c = coords.get(r["GeneSymbol"])
            if not c or r["SNPChr"] != c["chrom"]:
                continue
            p = int(r["SNPPos"])
            if int(c["start"]) - FLANK <= p <= int(c["end"]) + FLANK:
                hits[r["GeneSymbol"]].append(r)
    return coords, hits


def _decompress_consumed(buf):
    """Decompress whole BGZF blocks from buf; return (data, bytes_consumed).

    bytes_consumed counts only complete blocks, so the next slice can resume exactly
    at a block boundary and no record is split across two fetches.
    """
    out = bytearray()
    pos = 0
    while pos + 18 <= len(buf):
        ok, bsize = T.is_bgzf(buf[pos:pos + 18])
        if not ok or pos + bsize > len(buf):
            break
        out += zlib.decompress(buf[pos:pos + bsize], 16 + zlib.MAX_WBITS)
        pos += bsize
    return bytes(out), pos


class ChromCtx(object):
    """Per-chromosome cached tabix index + EUR column indices."""

    def __init__(self, chrom, keep):
        self.chrom = chrom
        self.url = KG % chrom
        self.sess = requests.Session()
        # the 1000genomes host intermittently closes connections; retry rather than
        # aborting a run that has already checkpointed earlier chromosomes
        last = None
        for attempt in range(5):
            try:
                ri = self.sess.get(self.url + ".tbi", timeout=600)
                ri.raise_for_status()
                break
            except Exception as exc:
                last = exc
                time.sleep(3 * (attempt + 1))
        else:
            raise RuntimeError("index fetch failed for chr%s: %r" % (chrom, last))
        BYTES["tbi"] += len(ri.content)
        self.tbi = T.Tbi(ri.content)
        for attempt in range(5):
            try:
                lines, st = T.header_lines(self.url, sess=self.sess, nbytes=1 << 18)
                break
            except Exception as exc:
                last = exc
                time.sleep(3 * (attempt + 1))
        else:
            raise RuntimeError("header fetch failed for chr%s: %r" % (chrom, last))
        BYTES["header"] += st["bytes"]
        chrom_line = next(l for l in lines if l.startswith("#CHROM"))
        samples = chrom_line.rstrip("\r").split("\t")[9:]
        self.n_samples_total = len(samples)
        self.idx = np.array([i for i, s in enumerate(samples) if s in keep])
        self.n_eur = len(self.idx)

    def region(self, start, end):
        """Fetch each index chunk as its own byte range.

        Collapsing disjoint chunks into a single min->max span is catastrophic here:
        1000G phase 3 VCFs carry structural variants that scatter a cis-window across
        15-17 disjoint chunks, and the collapsed span reached 491 MB for a 230 kb window
        (vs 2.5 MB when the window happens to land in one chunk). Results were unaffected
        -- the extra records get filtered out by position -- but the fetch never finished.
        """
        ch = self.tbi.chunks(self.chrom, start - 1, end)
        if not ch:
            return []
        GAP = 1 << 16
        spans = []
        for cb, ce in sorted(ch):
            lo, hi = cb >> 16, (ce >> 16) + T.MAX_BGZF_BLOCK
            if spans and lo - spans[-1][1] <= GAP:
                spans[-1][1] = max(spans[-1][1], hi)
            else:
                spans.append([lo, hi])

        SLICE = 1 << 20
        out = []
        for lo, hi in spans:
            cur = lo
            carry = b""
            while cur < hi:
                stop = min(cur + SLICE - 1, hi)
                stats = {}
                buf, _ = T.http_get_range(self.url, cur, stop, self.sess, stats)
                BYTES["data"] += stats.get("bytes", 0)
                if not buf:
                    break
                data, consumed = _decompress_consumed(buf)
                if consumed == 0:
                    break
                # A VCF record can span a BGZF block boundary, so resuming at the next
                # block splits that record and silently loses it. Carry the trailing
                # partial line into the next slice; consecutive BGZF blocks concatenate
                # back into one continuous byte stream, so this reconstructs it exactly.
                data = carry + data
                lines = data.split(b"\n")
                carry = lines.pop()
                out.extend(lines)
                cur += consumed
                # records inside a chunk are position-sorted, so once we see one past
                # the region end this chunk has nothing further to contribute
                past = False
                for l in reversed(lines):
                    if not l or l[:1] == b"#":
                        continue
                    p = l.split(b"\t", 3)
                    if len(p) < 3:
                        continue
                    try:
                        past = int(p[1]) > end
                    except ValueError:
                        continue
                    break
                if past:
                    break
        return out


def dosages(ctx, start, end, want):
    """want: dict (pos, frozenset(alleles)) -> key. Returns {key: np.array}."""
    out = {}
    for raw in ctx.region(start, end):
        if not raw or raw[:1] == b"#":
            continue
        f = raw.split(b"\t", 10)
        if len(f) < 10:
            continue
        try:
            pos = int(f[1])
        except ValueError:
            continue
        if not (start <= pos <= end):
            continue
        ref = f[3].decode()
        alt = f[4].decode()
        k = (pos, frozenset((ref, alt)))
        if k not in want:
            continue
        gts = raw.split(b"\t")[9:]
        # A byte-range fetch ends at a BGZF block boundary, which can fall mid-record,
        # so every span can yield one truncated line. Such a line may still parse a
        # position and allele pair while carrying too few genotype columns. Drop any
        # record whose genotype count does not match the VCF header exactly.
        if len(gts) != ctx.n_samples_total:
            continue
        blob = b"".join(g[:3] for g in gts)
        if len(blob) != 3 * len(gts):
            continue
        arr = np.frombuffer(blob, dtype="S1").reshape(-1, 3)
        d = (arr[:, 0] == b"1").astype(np.float32) + (arr[:, 2] == b"1").astype(np.float32)
        out[want[k]] = d[ctx.idx]
    return out


def clump_one(ctx, sym, coord, rows):
    rows = sorted(rows, key=lambda r: float(r["Pvalue"]))
    start = int(coord["start"]) - FLANK
    end = int(coord["end"]) + FLANK
    want = {}
    for r in rows:
        want[(int(r["SNPPos"]), frozenset((r["AssessedAllele"], r["OtherAllele"])))] = r["SNP"]
    dos = dosages(ctx, start, end, want)

    keys, vecs = [], []
    for r in rows:
        v = dos.get(r["SNP"])
        if v is None:
            continue
        keys.append(r)
        vecs.append(v)
    n_matched = len(keys)
    if not keys:
        return dict(symbol=sym, n_hits=len(rows), n_in_ld_panel=0, n_independent=0,
                    best_p=float(rows[0]["Pvalue"]), best_snp=rows[0]["SNP"],
                    independent_snps="", note="no_hit_snp_found_in_1000G_panel")

    X = np.vstack(vecs)
    sd = X.std(axis=1)
    ok = sd > 0
    mu = X.mean(axis=1, keepdims=True)
    Z = np.zeros_like(X)
    Z[ok] = (X[ok] - mu[ok]) / sd[ok][:, None]
    n = X.shape[1]

    kept = []
    for i in range(len(keys)):
        if not ok[i]:
            continue
        if kept:
            r = (Z[kept] @ Z[i]) / n
            if np.max(r * r) >= R2_CLUMP:
                continue
        kept.append(i)
    return dict(symbol=sym, n_hits=len(rows), n_in_ld_panel=n_matched,
                n_independent=len(kept), best_p=float(keys[0]["Pvalue"]),
                best_snp=keys[0]["SNP"],
                independent_snps=";".join(keys[i]["SNP"] for i in kept), note="")


def main():
    keep = eur_samples()
    coords, hits = load_inputs()
    print("EUR samples in 1000G phase3 panel:", len(keep))
    print("targets to clump:", len(hits))

    by_chrom = collections.defaultdict(list)
    for sym in hits:
        by_chrom[coords[sym]["chrom"]].append(sym)

    os.makedirs(CKPT, exist_ok=True)
    cols = ["symbol", "n_hits", "n_in_ld_panel", "n_independent", "best_p", "best_snp",
            "independent_snps", "note"]
    # PER-TARGET checkpointing. Per-chromosome was too coarse: a chromosome that takes
    # longer than the run window never checkpointed, so every restart redid it from
    # scratch and made zero net progress.
    results = []
    t0 = time.time()
    lock = threading.Lock()
    todo = sorted(by_chrom, key=lambda c: (len(c), c))
    for chrom in todo:
        cpath = os.path.join(CKPT, "chr%s.tsv" % chrom)
        done = {}
        if os.path.exists(cpath):
            with open(cpath, encoding="utf-8") as fh:
                for r in csv.DictReader(fh, delimiter=chr(9)):
                    done[r["symbol"]] = r
        syms = [s for s in by_chrom[chrom] if s not in done]
        results.extend(done.values())
        if not syms:
            print("  chr%-3s complete from checkpoint (%d targets)" % (chrom, len(done)),
                  flush=True)
            continue
        fresh = not os.path.exists(cpath)
        fh = open(cpath, "a", newline="", encoding="utf-8")
        w = csv.DictWriter(fh, fieldnames=cols, delimiter=chr(9), lineterminator=chr(10))
        if fresh:
            w.writeheader()
            fh.flush()
        ctx = ChromCtx(chrom, keep)
        n_new = 0
        try:
            with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
                futs = {ex.submit(clump_one, ctx, s, coords[s], hits[s]): s for s in syms}
                for fu in cf.as_completed(futs):
                    r = fu.result()
                    with lock:
                        w.writerow(r)
                        fh.flush()
                        os.fsync(fh.fileno())
                        results.append(r)
                        n_new += 1
        finally:
            fh.close()
        print("  chr%-3s %4d new (+%d cached) n_eur=%d cum=%d %.0fs"
              % (chrom, n_new, len(done), ctx.n_eur, len(results), time.time() - t0),
              flush=True)

    results.sort(key=lambda r: r["symbol"])
    with open(D + r"\clumped_instruments.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(results)
    json.dump(dict(BYTES), open(D + r"\clump_bytes.json", "w"), indent=2)
    print("bytes over wire:", dict(BYTES), "total=%d" % sum(BYTES.values()))
    print("seconds: %.0f" % (time.time() - t0))
    print("wrote", D + r"\clumped_instruments.tsv")


if __name__ == "__main__":
    main()
