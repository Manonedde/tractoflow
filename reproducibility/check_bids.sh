#!/bin/bash
# Check that Tractoflow and sf-tractomics read the same inputs from a BIDS
# dataset. Light (no computation): run on a login node.
#   ./check_bids.sh <bids_dir> <work_dir>
# Prints one row per session: what each pipeline selected (DWI, reverse
# b0/DWI, T1, readout) and whether topup will run.
set -euo pipefail
BIDS=$1; W=$2
HERE=$(cd "$(dirname "$0")" && pwd)
IMG=${IMG_DIR:?}/scilus_2.2.2.sif
SFTRACTOMICS_DIR=${SFTRACTOMICS_DIR:?}
module load apptainer nextflow
export NXF_APPTAINER_CACHEDIR=$IMG_DIR MPLCONFIGDIR=$W
mkdir -p "$W"

# Tractoflow: Read_BIDS runs scil_bids_validate
apptainer exec --no-home -B /project -B /scratch -B /home "$IMG" \
    scil_bids_validate "$BIDS" "$W/tractoflow_bids_struct.json" --readout 0.062 -f -v WARNING \
    2> "$W/tractoflow_bids.log" || { cat "$W/tractoflow_bids.log"; exit 1; }

# sf-tractomics: stub run (BIDS parsing + process graph, files are touched only)
(cd "$W" && nextflow run "$SFTRACTOMICS_DIR/main.nf" -stub -profile full_reproducible,apptainer \
    --input "$BIDS" --outdir "$W/sf_stub" --run_gibbs_correction false \
    --run_t1_denoising false -with-trace "$W/sf_trace.txt" > "$W/sf_stub.log" 2>&1) \
    || { tail -30 "$W/sf_stub.log"; exit 1; }

python3 - "$W" <<'EOF'
import csv, glob, json, os, re, sys
W = sys.argv[1]
tf = {}
for e in json.load(open(f"{W}/tractoflow_bids_struct.json")):
    sid = f"sub-{e['subject']}" + (f"_ses-{e['session']}" if e.get("session") else "")
    rev = e.get("rev_dwi") or e.get("rev_topup") or ""
    tf[sid] = dict(dwi=os.path.basename(e.get("dwi") or ""), rev=os.path.basename(rev),
                   t1=os.path.basename(e.get("t1") or ""), readout=e.get("TotalReadoutTime"),
                   topup=bool(rev))
sf = {}
for t in csv.DictReader(open(f"{W}/sf_trace.txt"), delimiter="\t"):
    step = t["name"].split(" (")[0].split(":")[-1]
    if step not in ("PREPROC_TOPUP", "PREPROC_EDDY"):
        continue
    for wd in glob.glob(f"{W}/work/{t['hash']}*"):
        for f in os.listdir(wd):
            m = re.match(r"(sub-[^_]+(?:_ses-[^_]+)?)", f)
            if m:
                d = sf.setdefault(m.group(1), {"topup": False, "eddy": False})
                d["topup" if step == "PREPROC_TOPUP" else "eddy"] = True
# work files named after the subject only (no session): not a session
sf = {k: v for k, v in sf.items()
      if "_ses-" in k or not any(o.startswith(k + "_ses-") for o in sf)}
print(f"{'session':28s} {'TF topup':9s} {'SF topup':9s} {'TF readout':10s}  TF dwi | rev | t1")
bad = 0
for sid in sorted(set(tf) | set(sf)):
    a, b = tf.get(sid), sf.get(sid)
    ok = a is not None and b is not None and a["topup"] == b["topup"]
    bad += not ok
    print(f"{sid:28s} {str(a and a['topup']):9s} {str(b and b['topup']):9s} "
          f"{str(a and a['readout']):10s}  "
          + (f"{a['dwi']} | {a['rev'] or '-'} | {a['t1']}" if a else "NOT READ BY TRACTOFLOW")
          + ("" if b else "   NOT READ BY SF-TRACTOMICS") + ("" if ok else "   <-- MISMATCH"))
print(f"\n{len(tf)} sessions Tractoflow, {len(sf)} sf-tractomics, {bad} mismatch(es)")
print(f"Logs: {W}/tractoflow_bids.log, {W}/sf_stub/pipeline_info/BIDS_logs.txt")
EOF
