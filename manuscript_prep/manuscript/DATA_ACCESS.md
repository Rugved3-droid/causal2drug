# DATA_ACCESS — recovery-map resources

Checked 2026-10-10, after `RECOVERY_RULE.md` was written and before any classification.
**Result on 2026-10-10: three of four resources obtained; the Fenland table could not be
obtained, and work stopped there, as the rule requires.**

**Update 2026-10-11: the project owner downloaded the Fenland workbook in a browser. All four
resources are now readable. Classification was run only after this point
(`recovery_scripts/build_recovery_map.py`, output `recovery_map_summary.txt`).**

| resource | source tried | result | file |
|---|---|---|---|
| **GTEx v8** significant cis-eQTL pairs, all tissues | `https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar` | **Downloaded and readable.** 1,562,828,800 bytes; MD5 of the downloaded file `01syFSvbITFrJQnEawr5mA==` equals the MD5 the server publishes; the archive lists 49 `signif_variant_gene_pairs` files and 49 `egenes` files, one pair per tissue. Only the file listing has been read | `recovery_data/GTEx_Analysis_v8_eQTL.tar` |
| **UKB-PPP**, Sun et al. 2023, supplementary tables | Europe PMC supplementary-files service for PMC10567551 | **Downloaded and readable.** `41586_2023_6592_MOESM3_ESM.xlsx`, 34,411,394 bytes, 31 sheets | `recovery_data/41586_2023_6592_MOESM3_ESM.xlsx` |
| **deCODE**, Ferkingstad et al. 2021, supplementary tables | `static-content.springer.com … 41588_2021_978_MOESM4_ESM.xlsx` | **Downloaded and readable.** 44,366,108 bytes, 19 sheets | `recovery_data/ferkingstad2021_MOESM4.xlsx` |
| **Fenland**, Pietzner et al. 2021, supplementary tables | PubMed Central author-manuscript file `NIHMS1866945-supplement-Supplementary_tables.xlsx`, downloaded by the project owner in a browser to the root of drive D and copied unchanged | **Obtained 2026-10-11 and readable.** 4,076,443 bytes, MD5 `ba9ab9207ecd551b0c21aea1bff702f8` (identical for the original and the copy). ST2: 10,674 pQTL rows with `cis.trans`, `Chr`, `Position`, `P.value`, Ensembl ID and symbol. ST3: 4,030 aptamers with gene chromosome, start and end. Build not stated in the workbook; ST3 coordinates are identical to Ensembl GRCh37.87 for 3,917 of 3,933 single-gene aptamers (`recovery_map_summary.txt` [R08]) | `recovery_data/pietzner2021_supplementary_tables.xlsx` |

## Fenland: what was tried on 2026-10-10 (scripted access, superseded by the update above)

| route | outcome |
|---|---|
| science.org supplementary zip (`science.abj1541_tables_s1_to_s15.zip`) | HTTP 403, with and without a browser user-agent |
| Europe PMC supplementary-files service for PMC9904207 | "Article with id PMC9904207 is not open access one" |
| PubMed Central author-manuscript file `NIHMS1866945-supplement-Supplementary_tables.xlsx` | The file exists and is linked from the PMC page, but PMC answers a scripted request with a "Preparing to download" page that requires a browser to solve a proof-of-work challenge |
| Europe PMC mirror of the same file | HTTP 403 |

The PMC challenge is an access control aimed at automated downloads. I did not try to get
around it. The file is free to download in an ordinary browser:

`https://pmc.ncbi.nlm.nih.gov/articles/instance/9904207/bin/NIHMS1866945-supplement-Supplementary_tables.xlsx`

**To unblock:** download that file in your browser and save it as
`D:\causal2drug\manuscript_prep\manuscript\recovery_data\pietzner2021_supplementary_tables.xlsx`.
Or tell me to proceed without Fenland, in which case every Fenland cell is recorded as
UNVERIFIED and "locked" is stated relative to GTEx v8, UKB-PPP and deCODE only.

I did not substitute another source for Fenland (for example the study's data portal or its
preprint).

## What the two readable pQTL workbooks contain

Only sheet names and header rows were read to establish this. Two data rows mentioning PCSK9
appeared in the header printout (a coding-variant table and an MR example in the UKB-PPP
workbook); neither is a cis-pQTL result and neither was used.

| study | sheet for the rule | columns that matter | consequence for the rule |
|---|---|---|---|
| UKB-PPP | ST9 (discovery) and ST10 (combined): significant pQTLs, 30,747 rows each. ST3: the 2,941 assayed proteins with gene symbol, UniProt, gene chromosome, start and end | `cis/trans`, `GENPOS (hg38)`, `log10(p)`, `Ensembl gene ID`, `Distance to gene` | Distance from the gene body can be checked with the study's own gene coordinates. "Assayed" can be read from ST3. Note the study's own threshold is p < 1.7 × 10⁻¹¹, stricter than 5 × 10⁻⁸, so variants between the two are not in the table |
| deCODE | ST02: sentinel pQTLs, 28,193 rows. ST01: the 5,284 aptamers with gene, UniProt and Ensembl ID | `cis/trans`, `chr (var.)`, `pos (var.)` (hg38), `TSS (prot.)`, `distTSS`, p-value columns | The table gives the TSS, not the gene end. The gene body will be taken from GENCODE v26 (GRCh38), the same build, as carried in the GTEx files. This is a clarification of the rule, not a change, and will be logged under Amendments when classification starts |

## Open point on the rule, for your decision before I classify

UKB-PPP reports pQTLs only below p < 1.7 × 10⁻¹¹ and deCODE below its own study-wide
threshold. A cis variant with 5 × 10⁻⁸ > p > the study threshold is therefore invisible in the
published tables. The rule already says the pQTL check is conservative for this reason. No
change is proposed; it is recorded so the limit is explicit.
