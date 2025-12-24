import pandas as pd, subprocess, json, os
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

MANIFEST="manifest_af3only.csv"
OUT="mmgbsa_openmm_scores_af3only.csv"
NWORKERS=60  # use your 64-core box

def run_one(path):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    single_script = os.path.join(script_dir, "mmgbsa_openmm_single.py")
    python_exe = "/home/kiran/miniforge3/envs/mmgbsa_openmm/bin/python"
    cmd=[python_exe, single_script, path]
    p=subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode!=0:
        return {"file": path, "status":"ERROR", "error": p.stderr[:200]}
    try:
        return json.loads(p.stdout.strip())
    except Exception:
        return {"file": path, "status":"ERROR", "error":"JSON_PARSE", "raw": p.stdout[:200]}

def main():
    df=pd.read_csv(MANIFEST)
    files=[x for x in df["file"].tolist() if isinstance(x,str) and x.endswith(".cif")]
    rows=[]
    with ProcessPoolExecutor(max_workers=NWORKERS) as ex:
        futs={ex.submit(run_one,f):f for f in files}
        for fut in tqdm(as_completed(futs), total=len(futs)):
            rows.append(fut.result())
    out=pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print("Wrote", OUT, "rows=", len(out))
    print(out["status"].value_counts(dropna=False))
if __name__=="__main__":
    main()
