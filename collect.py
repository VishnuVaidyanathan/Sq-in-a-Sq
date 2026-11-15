"""Gather search results: keep the smallest checked box for each n and write results/summary.csv.

    python collect.py incoming_folder
"""
import csv, glob, math, os, sys
import numpy as np
import band_table as bt
import search


def checked(path):
    v, side, n = search.load(path)
    T = np.c_[v[:n], v[n:2*n], np.degrees(v[2*n:])]
    return n, side, bt.check(T, side), v


def main(incoming):
    os.makedirs("results", exist_ok=True)
    kept = {}
    for path in sorted(glob.glob("results/*.csv")) + sorted(glob.glob(os.path.join(incoming, "**", "*.csv"), recursive=True)):
        if os.path.basename(path) == "summary.csv": continue
        try:
            n, side, overlap, v = checked(path)
        except Exception as e:
            print("skipped", path, e); continue
        if overlap > 1e-9:
            print("rejected (tiles overlap)", path, overlap); continue
        if n not in kept or side < kept[n][0] - 1e-9: kept[n] = (side, v, open(path).read())
    if not kept: print("no results"); return
    best = bt.build_best(max(kept))
    with open("results/summary.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["n", "as m^2+x", "band or grid side", "search side", "gain"])
        for n in sorted(kept):
            side, v, text = kept[n]; m = math.isqrt(n); base = best[n][0]
            if side < base - 1e-7:
                open(f"results/{n}.csv", "w").write(text)
            elif os.path.exists(f"results/{n}.csv"):
                os.remove(f"results/{n}.csv")          # a result that does not beat the band table is not kept
            w.writerow([n, f"{m}^2+{n - m*m}", f"{base:.6f}", f"{side:.6f}", f"{base - side:.6f}"])
            print(n, f"band/grid {base:.6f}  search {side:.6f}  gain {base - side:.6f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "incoming")
