"""AUDIT extra checks, not rules of the project.

[B09] status by highest phase (KEY_NUMBERS section 3).
[F01] are any targets recorded as 'not assayed' in a proteomic study present in that study's protein list
      under an HGNC previous symbol or a UniProt accession? (RECOVERY_RULE clarification 5 matches on the
      table's symbol or Ensembl ID only.)
Run:  python audit/a08_extra_checks.py
"""
import csv, collections

O = r"D:\causal2drug\audit\out" + "\\"


def rd(f, skip):
    r = list(csv.reader(open(f, encoding="utf-8"), delimiter="\t"))
    return r[skip], r[skip + 1:]


r = list(csv.DictReader(open(O + "audit_per_target.csv")))
for ph, name in zip("12345", ["Preclinical", "Phase I", "Phase II", "Phase III", "Launched"]):
    x = [t for t in r if t["phase"] == ph]
    c = collections.Counter(t["status"] for t in x); s = collections.Counter(t["subtype"] for t in x)
    na = sum(1 for t in x if t["chrom"] not in ("X", "Y", "MT"))
    print(f"[B09] {name}: {len(x)} targets; instrumented {c['INSTRUMENTED']}; no eQTL data {c['NO_EQTL_DATA']} "
          f"({s['autosomal']} / {s['x_linked']} + {s['mitochondrial']}); autosomal-only {s['autosomal']} of {na}; "
          f"no qualifying instrument {c['NO_QUALIFYING_INSTRUMENT']}; upper {c['NO_EQTL_DATA']+c['NO_QUALIFYING_INSTRUMENT']}; "
          f"unresolved {c['UNRESOLVED']}")

P = O + "pqtl\\"
hg = {x["symbol"]: x for x in csv.DictReader(open(O + "hgnc_complete_set.txt", encoding="utf-8"), delimiter="\t")}
h, d = rd(P + "ukbppp_ST3.tsv", 2)
ukb_s = {s.strip() for x in d for s in x[4].split(";")}; ukb_u = {u.strip() for x in d for u in x[5].replace(";", "|").split("|")}
h, d = rd(P + "decode_ST01.tsv", 2)
dec_s = {s.strip() for x in d for s in x[3].split("|")}; dec_u = {u.strip() for x in d for u in x[6].split("|")}
h, d = rd(P + "fenland_ST3.tsv", 0)
fen_s = {s.strip().lstrip("'") for x in d for s in x[16].split("|")}; fen_u = {u.strip() for x in d for u in x[13].replace(",", "|").split("|")}
hits = []
for t in csv.DictReader(open(O + "audit_recovery_map.csv")):
    g = hg.get(t["symbol"])
    if not g:
        continue
    prev = set(g["prev_symbol"].strip('"').split("|")) - {""}; uni = set(g["uniprot_ids"].strip('"').split("|")) - {""}
    for study, S, U, flag in (("UKB-PPP", ukb_s, ukb_u, t["ukb_assayed"]), ("deCODE", dec_s, dec_u, t["decode_listed"]),
                              ("Fenland", fen_s, fen_u, t["fenland_assayed"])):
        if flag == "False" and ((prev & S) or (uni & U)):
            hits.append((t["symbol"], t["tier"], study, sorted(prev & S), sorted(uni & U)))
print("[F01] recorded not assayed, but present under a previous symbol or UniProt accession:", len(hits))
for x in hits:
    print("      ", x)
print("[F02] of these, currently 'locked: never assayed':", [x for x in hits if x[1] == "locked: never assayed"])
