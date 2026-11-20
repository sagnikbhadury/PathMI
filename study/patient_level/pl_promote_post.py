"""pl_promote_post.py -- after pl_promote.R:
1. converts the CSV sources into the three merged workbooks in TRIDENT/ (openpyxl);
2. collapses each submeans.npy (n_slides x 16 x 20 x D, float32) to patients
   (mean over a patient's slides of the m-tile means), written to the new pool_pass.
"""
import csv, os, sys, numpy as np, openpyxl
TRI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(TRI, "patient_level", "patient", "xlsx_src")
for name in ("mi summary table", "best kpca summary", "incremental mi summary"):
    wb = openpyxl.Workbook(); ws = wb.active
    with open(os.path.join(SRC, name + ".csv"), newline="", encoding="utf-8") as f:
        for i, row in enumerate(csv.reader(f)):
            if i == 0: ws.append(row); continue
            out = []
            for v in row:
                try: out.append(float(v) if v not in ("", "NA") else None)
                except ValueError: out.append(v)
            ws.append(out)
    wb.save(os.path.join(TRI, name + ".xlsx")); print("wrote", name + ".xlsx")

POOL_SLIDE = "M:/bhadury/mi_results/pool_pass_slide"; POOL = "M:/bhadury/mi_results/pool_pass"
import subprocess
RS = "C:/Users/bhadury/AppData/Local/Programs/R/R-4.4.3/bin/Rscript.exe"
for tag in sorted(os.listdir(POOL_SLIDE)):
    src = os.path.join(POOL_SLIDE, tag, "submeans.npy"); dst = os.path.join(POOL, tag, "submeans.npy")
    if not os.path.exists(src) or os.path.exists(dst): continue
    # slide roster order from the slide pooled.rds via a tiny R call
    r = subprocess.run([RS, "-e", "cat(readRDS('%s/%s/pooled.rds')$roster, sep='\\n')" % (POOL_SLIDE, tag)], capture_output=True, text=True)
    roster = [s for s in r.stdout.split("\n") if s]
    pat = [s[:12] for s in roster]; units = list(dict.fromkeys(pat))
    S = np.load(src, mmap_mode="r"); assert S.shape[0] == len(roster), (tag, S.shape, len(roster))
    out = np.zeros((len(units),) + S.shape[1:], dtype=np.float32)
    idx = {u: [i for i, p in enumerate(pat) if p == u] for u in units}
    for j, u in enumerate(units): out[j] = np.asarray(S[idx[u]], dtype=np.float64).mean(axis=0)
    np.save(dst, out); print(tag, S.shape, "->", out.shape)
