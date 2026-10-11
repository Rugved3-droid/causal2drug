# Deliverable 0 — Access verification

Working dir: `D:\causal2drug` · All commands run 2026-08-26 UTC · Python 3.12.10, Git Bash on Windows 11

**Toolchain constraint (verified, not assumed).** `tabix` and `bgzip` binaries are absent and
`pysam` has no wheel for this platform:

```bash
for t in curl tabix bgzip python; do printf "%-8s " "$t"; command -v "$t" || echo MISSING; done
python -m pip install --only-binary=:all: pysam    # -> ERROR: No matching distribution found for pysam
```

So region queries are served by a from-scratch BGZF + tabix-index + HTTP-Range client:
`src/tabix_http.py`. It parses the `.tbi` binary format, computes candidate bins via the htslib
`reg2bins` scheme, bounds the byte span with the linear index, issues one ranged GET, and reports
exactly how many bytes crossed the wire. Verifying with this client is a stronger test than
verifying with `tabix` itself, because it proves the *remote* file is range-servable rather than
trusting a local wrapper.

---

## 1. GWAS Catalog harmonised summary statistics — VERIFIED

Target: `GCST90239658` (Graham SE 2021, LDL cholesterol, max N = 1,320,016), selected mechanically
(see §"How the study was chosen"). Data file is **1,416,980,601 bytes**; sibling `.tbi` is
**1,803,528 bytes**.

Region: HMGCR ± 100 kb = `chr5:75,236,329-75,464,241`. Gene bounds from Ensembl REST (GRCh38):
`ENSG00000113161 chr5:75,336,329-75,364,241`.

```bash
cd /d/causal2drug
python src/tabix_http.py \
  "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90239001-GCST90240000/GCST90239658/harmonised/GCST90239658.h.tsv.gz" \
  5 75236329 75464241
```

| quantity | value |
|---|---|
| `Accept-Ranges` on data file | `bytes` |
| HTTP status of ranged GETs | `206 Partial Content` |
| data file is BGZF (not plain gzip) | `true`, first block 24,514 B |
| byte span requested | `[472,443,986 – 472,603,180]` |
| **bytes transferred total** | **1,962,741** (of which 1,803,528 is the one-off `.tbi`) |
| **fraction of the 1.42 GB file pulled** | **0.1385 %** |
| variants returned in window | **3,323** |
| wall clock | 5.14 s |

All 3,323 returned positions lie inside the requested window and are ascending
(`pos_min=75,236,454 pos_max=75,464,115`). Extract + provenance saved to
`deliverables/d0/extracts/GWASCAT_GCST90239658_chr5_75236329_75464241.tsv`.

**Caveat that matters — `.tbi` is NOT universal.** The index is a per-study artefact, not a
platform guarantee. Older studies lack it entirely; e.g. `GCST002222` (Willer 2013 GLGC LDL) has
only `.h.tsv.gz`, `-meta.yaml`, `readme.txt` and no `.tbi`, so it is **not** range-queryable.

Measured by HEAD-probing a seeded random sample of the 156,116-line harmonised manifest:

```bash
python src/check_tbi.py --list data/raw/harmonised_list.txt \
  --out deliverables/d0/tbi_probe_random200.csv --n 200 --seed 20260826
```

* random sample n = 200: **173 have `.tbi` → 86.50 %**, Wald 95 % CI **[81.8 %, 91.2 %]**
* LDL-trait subset, **all 294 probed** (census, no sampling error): **219 → 74.49 %**
* data-file sizes in the random sample: min 2.72 MB, median **253 MB**, max **4.17 GB**

## 2. eQTL Catalogue — VERIFIED

Release R8, dataset `QTD000266` = GTEx_v10 liver, sample_size 269, `quant_method=ge`
(from `data_tables/dataset_metadata_r8.tsv`). Data file **2,906,926,368 bytes**, `.tbi`
**2,112,225 bytes**. Every dataset directory carries `.all.tsv.gz{,.tbi}` and `.cc.tsv.gz{,.tbi}`.

```bash
python src/tabix_http.py \
  "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000266/QTD000266.all.tsv.gz" \
  5 75236329 75464241
```

| quantity | value |
|---|---|
| `Accept-Ranges` / status | `bytes` / `206` |
| BGZF | `true`, first block 11,126 B |
| byte span | `[2,219,071,243 – 2,219,315,943]` |
| **bytes transferred** | **2,356,944** (2,112,225 of it the `.tbi`) |
| **fraction of the 2.91 GB file** | **0.0811 %** |
| rows in window (all molecular traits) | 8,990 across 16 distinct traits |
| rows for HMGCR `ENSG00000113161` | **592** |
| wall clock | 6.02 s |

Extract saved to `.../EQTLCAT_QTD000266_chr5_75236329_75464241_HMGCR.tsv`.

**The eQTL Catalogue REST API is retired** — `https://www.ebi.ac.uk/eqtl/api/v2/datasets` returns
`HTTP 410`: *"This API is no longer available."* The FTP + tabix path verified above is the
supported route.

## 3. Sources requiring full-file download

All sizes below are measured HEAD `Content-Length` values, not estimates.

**Do NOT require full download (range-queryable, verified above):**

| source | full size | pulled per cis-window query |
|---|---|---|
| GWAS Catalog harmonised sumstats (when `.tbi` present, ~86.5 %) | 253 MB median / 4.17 GB max | ~160 KB + 1.8 MB index |
| eQTL Catalogue `.all.tsv.gz` | 2.91 GB (this dataset) | ~245 KB + 2.1 MB index |
| 1000 Genomes GRCh38 phased VCFs (LD panel) | 30.20 GiB whole release | ~2.24 MB + 179 KB index |

The LD reference panel being range-queryable was not assumed — it was tested, and the test failed
first (see §5), which is why the number is trustworthy now:

```bash
python src/tabix_http.py \
  ".../20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr5.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" \
  chr5 75236329 75464241
# 5,644 records, 3,202 samples, 2,422,986 bytes = 0.1276 % of the 1.90 GB chr5 file
```

**DO require full download:**

| source | measured bytes | note |
|---|---|---|
| Minikel `pp.tsv` | 3,162,501 | downloaded, verified — §5 |
| GWAS Catalog `harmonised_list.txt` | 11,959,625 | downloaded |
| GWAS Catalog studies release v1.0.3.1 | 129,139,393 | downloaded; needed for trait→study join |
| eQTL Catalogue `dataset_metadata_r8.tsv` | 123,426 | downloaded |
| NLM MeSH `desc2026.xml` | 312,952,703 | not yet downloaded; `Accept-Ranges: bytes` |
| EFO `efo.obo` (latest release) | 88,382,139 | not yet downloaded |
| MONDO `mondo.obo` | 53,086,111 | not yet downloaded |
| Open Targets `target/` (parquet) | 89,931,558 | 3 files |
| Open Targets `target_prioritisation/` | 1,915,147 | tractability/safety scores |
| Open Targets `disease/` | 7,150,787 | |
| Open Targets `drug_molecule/` | 14,647,216 | |
| Open Targets `drug_mechanism_of_action/` | 581,210 | |
| Open Targets `drug_warning/` | 255,356 | |

Sum of the not-yet-downloaded items ≈ **568 MB**. Current project footprint: **141 MB**;
free space on D: **846 GB**.

**Open Targets release caveat.** `latest` is byte-identical to release `26.06` (same parquet part
UUIDs). The release-level `manifest.json` for `26.06` reports `"result": "failure"`
(`started_at 2026-06-15`); `25.12` also reports `failure` and `26.03` reports `pending`. Individual
dataset directories do carry `_SUCCESS` markers. I have **not** verified what the release-level
failure covers. Recommendation: pin an explicit release rather than `latest`, and check per-dataset
`_SUCCESS` before use. Flagged, not resolved.

## 4. Credentials

**Nothing on the verified critical path needs a credential.** GWAS Catalog FTP + REST, eQTL
Catalogue FTP, 1000 Genomes FTP, Open Targets FTP, EBI OLS4, Ensembl REST and the Minikel GitHub
repo all answered `200` unauthenticated.

Two sources **do** demand credentials. Neither is currently required, so I am **not** working
around them and **not** asking you to supply them yet:

| source | status | would be needed if |
|---|---|---|
| **OpenGWAS** `api.opengwas.io` | `HTTP 401` — *"From 1st May 2024 you must provide a token (JWT)"* | we used OpenGWAS for GWAS retrieval or its LD-clumping endpoint instead of the verified GWAS Catalog + local 1000G route |
| **UMLS UTS** `uts-ws.nlm.nih.gov` | `HTTP 401` — API key required | Deliverable 3's MeSH↔EFO mapping cannot be done from open MeSH/EFO/MONDO files alone |

Synapse: **not required** by any source in the verified stack.

**One thing to decide before Deliverable 3**, so flagging now rather than at the point of need:
`pp.tsv` keys indications by **MeSH descriptor ID** (`D002545` etc.), while GWAS Catalog keys traits
by **EFO**. Whether that crosswalk is achievable from open files only (MeSH XML + EFO + MONDO
xrefs) or needs a UMLS key is a Deliverable 3 question. I have not tested it, so it is **UNVERIFIED**.

## 5. `pp.tsv` checksum and contents

**Checksum — VERIFIED, with an important qualification about what "published release" means here.**

The repo `ericminikel/genetic_support` is CC-BY-4.0 (confirmed via GitHub API), 443,617 KB, and has
exactly **one** release: `v1.0`, published `2024-03-05T18:18:26Z`. That release ships **no assets**
and the repository contains **no `md5sum`/`sha256sum`/checksum manifest** — verified by scanning all
192 entries of the git tree. So there is no author-published checksum to compare against.

The one publisher-side integrity anchor that does exist is the git object hash, which GitHub
advertises independently of the file bytes. `pp.tsv` at tag `v1.0` and at `main` carry the **same**
blob SHA-1, so the file is unchanged since the release.

```bash
curl -s "https://api.github.com/repos/ericminikel/genetic_support/contents/data/pp.tsv?ref=v1.0" \
  | python -c "import sys,json;d=json.load(sys.stdin);print(d['size'],d['sha'])"
# 3162501 16d641186c391e188623c757a2feecef333d4a85

curl -sS -o pp.tsv "https://raw.githubusercontent.com/ericminikel/genetic_support/v1.0/data/pp.tsv"
python -c "
import hashlib; b=open('pp.tsv','rb').read()
print(hashlib.sha1(b'blob %d\x00'%len(b)+b).hexdigest())"
# 16d641186c391e188623c757a2feecef333d4a85  -> MATCH
```

| | |
|---|---|
| bytes | 3,162,501 (== advertised size) |
| git blob SHA-1 | `16d641186c391e188623c757a2feecef333d4a85` — **matches GitHub** |
| SHA-256 | `db1bda5367759fb223f9d6c4dcb8bf13eaf50229289c5f6ac152b4478f6e4fdb` |
| MD5 | `f994da2c7069528746e11602a2e0895e` |

The SHA-256 and MD5 are ours, recorded for future re-verification; there is nothing upstream to
compare them to. **UNVERIFIED:** whether the Nature 2024 paper or its supplement publishes an
independent checksum for `pp.tsv`. I did not access the paper.

**Shape.** 3,162,501 bytes, LF line endings, no CR bytes, trailing newline present,
29,477 lines = 1 header + **29,476 data rows**, all rows exactly 19 fields.
2,517 distinct genes · 1,014 distinct MeSH indications.

Columns: `ti_uid, gene, indication_mesh_id, indication_mesh_term, hcat, acat, ccat, hcatnum,
acatnum, ccatnum, succ_p_1, succ_1_2, succ_2_3, succ_3_a, oto_hcat, oto_acat, oto_ccat, orphan,
year_launch`

First 20 data rows: see `deliverables/d0/pp_first20.txt` (printed verbatim in the session report).

**Duplicates — 16 duplicated `ti_uid` keys, and they are not all benign.**

```bash
python src/pp_duplicates.py    # regenerates everything below from data/raw/pp.tsv
```

* exact duplicate **full lines**: 9 distinct lines, 9 excess rows
* duplicated **`ti_uid`**: **16** keys (29,460 distinct of 29,476 rows)
* duplicated **(gene, indication_mesh_id)**: 16 — identical set, confirming `ti_uid` is just
  `gene-mesh_id` concatenated

Of the 16, **4 carry conflicting outcome values** — the same target-indication program appears
twice with a *different highest phase reached*:

| ti_uid | conflict |
|---|---|
| `DAO-D022124` (Hyperammonemia) | Phase III vs Preclinical |
| `HRH1-D000086382` (COVID-19) | **Launched** vs Phase II |
| `IFNAR2-D012141` (Respiratory Tract Infections) | **Launched** vs Phase II |
| `MC1R-D010787` (Photosensitivity Disorders) | Phase III vs Phase II |

A further 3 (`MAPK14-`, `TGFB1-`, `TNFRSF1A-D000086382`) differ only in `oto_*`/`orphan` columns,
and the remaining 9 are byte-identical.

Dedup consequences: **29,467** rows if deduplicating identical lines only; **29,460** if
deduplicating by `ti_uid`. The 4 outcome-conflicting keys need an explicit, pre-registered
resolution rule — and because two of them flip between *Launched* and *Phase II*, the choice is
outcome-bearing. It must be fixed before any outcome is looked at. Not decided here.

---

### How the GWAS study was chosen (mechanical, reproducible)

```bash
# accessions with harmonised sumstats
sed 's|^\./||' data/raw/harmonised_list.txt \
  | sed 's|.*/\(GCST[0-9]*\)/harmonised/.*|\1|' | sort -u   # 155,552 accessions

# studies mapped to EFO_0004611 (LDL cholesterol measurement, label from GWAS Catalog API)
awk -F'\t' 'NR>1 && $14 ~ /EFO_0004611/ {print $15"\t"$3"\t"$4"\t"$9}' \
  data/raw/gwas-catalog-download-studies-v1.0.3.1.txt | sort -u   # 740 studies

# intersect -> 292 with harmonised sumstats -> probe for .tbi -> rank by sample size
```

`GCST90239658` is the largest such study that has a `.tbi`.

### Files produced

```
D:\causal2drug\
  src\tabix_http.py         pure-Python BGZF+tabix HTTP-range client
  src\check_tbi.py          .tbi availability prober (seeded sampling)
  src\dirsize.py            FTP directory size measurement
  src\pp_duplicates.py      pp.tsv duplicate analysis
  data\raw\                 5 files, read-only (chmod 444), 138 MB
  data\derived\             6 intermediate join files
  logs\retrieval_log.csv    url / timestamp / bytes / sha256 / md5 per raw file
  deliverables\d0\          this report, extracts, tbi probe CSVs, pp_first20.txt
```
