"""Build and check 45-degree band packings of unit squares in a square.

A layout = one rectangular block of p x q unit squares turned 45 degrees (long side along the
main diagonal), plus straight unit squares stacked from the four corners.
Every layout that is reported is checked tile against tile.
"""
import math, numpy as np
R2 = math.sqrt(2.0); EPS = 1e-9

def block_tiles(s, p, q, delta, t):
    cx = s/2 - delta/2 + t; cy = s/2 + delta/2 + t
    u = np.array([1, 1])/R2; v = np.array([-1, 1])/R2
    a = np.arange(q) - (q-1)/2; b = np.arange(p) - (p-1)/2
    A, B = np.meshgrid(a, b, indexing="ij")
    P = np.array([cx, cy]) + A.reshape(-1, 1)*u + B.reshape(-1, 1)*v
    return P, (cx, cy)

def corner_cells(s, corner, p, q, c):
    """Straight cells stacked from one corner that miss the block, kept as a staircase."""
    K = int(math.floor(s + EPS))
    i, j = np.meshgrid(np.arange(K), np.arange(K), indexing="ij")
    x0 = i.astype(float) if corner in ("BL", "TL") else s - i - 1.0
    y0 = j.astype(float) if corner in ("BL", "BR") else s - j - 1.0
    dx = x0 + .5 - c[0]; dy = y0 + .5 - c[1]
    Rb = (p + q)/(2*R2)
    hit = (np.abs(dx) < .5 + Rb - EPS) & (np.abs(dy) < .5 + Rb - EPS) & \
          (np.abs(dx + dy)/R2 < q/2 + R2/2 - EPS) & (np.abs(dy - dx)/R2 < p/2 + R2/2 - EPS)
    ok = ~hit
    ok = np.logical_and.accumulate(ok, axis=0); ok = np.logical_and.accumulate(ok, axis=1)
    return np.c_[x0[ok], y0[ok]]

def build(s, p, q, delta, t, order=("TL", "BR", "BL", "TR")):
    P, c = block_tiles(s, p, q, delta, t)
    # block must lie inside the box
    h = R2/2
    if P[:, 0].min() - h < -EPS or P[:, 1].min() - h < -EPS or P[:, 0].max() + h > s + EPS or P[:, 1].max() + h > s + EPS:
        return None
    acc = np.zeros((0, 2))
    for cn in order:
        cand = corner_cells(s, cn, p, q, c)
        if len(acc) and len(cand):
            clash = (np.abs(cand[:, None, 0] - acc[None, :, 0]) < 1 - EPS) & (np.abs(cand[:, None, 1] - acc[None, :, 1]) < 1 - EPS)
            cand = cand[~clash.any(axis=1)]
        acc = np.vstack([acc, cand])
    # tiles as (centre x, centre y, angle in degrees)
    T = np.vstack([np.c_[P, np.full(len(P), 45.0)], np.c_[acc + .5, np.zeros(len(acc))]])
    return T

def add_L(T, s):
    K = int(math.floor(s + EPS))
    top = np.c_[np.arange(K + 1) + .5, np.full(K + 1, s + .5), np.zeros(K + 1)]
    right = np.c_[np.full(K, s + .5), np.arange(K) + .5, np.zeros(K)]
    return np.vstack([T, top, right]), s + 1.0

def check(T, s):
    """Largest overlap between any two tiles, or of a tile with a wall. 0 means a valid packing."""
    x, y, a = T[:, 0], T[:, 1], np.radians(T[:, 2])
    h = .5*(np.abs(np.cos(a)) + np.abs(np.sin(a)))
    wall = max((h - x).max(), (x + h - s).max(), (h - y).max(), (y + h - s).max(), 0.0)
    worst = 0.0; n = len(T); B = 400
    for i0 in range(0, n, B):
        xi, yi, ai = x[i0:i0+B, None], y[i0:i0+B, None], a[i0:i0+B, None]
        dx = x[None, :] - xi; dy = y[None, :] - yi
        near = (dx*dx + dy*dy < 2.0001)
        near[np.arange(len(xi)), np.arange(i0, i0+len(xi))] = False
        I, J = np.nonzero(near)
        if not len(I): continue
        ddx, ddy = dx[I, J], dy[I, J]; a1 = a[i0 + I]; a2 = a[J]; d = a1 - a2
        hh = .5 + .5*(np.abs(np.cos(d)) + np.abs(np.sin(d)))
        o = np.minimum.reduce([hh - np.abs(ddx*np.cos(a1) + ddy*np.sin(a1)), hh - np.abs(-ddx*np.sin(a1) + ddy*np.cos(a1)),
                               hh - np.abs(ddx*np.cos(a2) + ddy*np.sin(a2)), hh - np.abs(-ddx*np.sin(a2) + ddy*np.cos(a2))])
        worst = max(worst, float(o.max()))
    return max(worst, float(wall))

def frac(z): return z - math.floor(z)

def search_block(m):
    """All band layouts with side between m and m+1. Returns list of (side, count, recipe)."""
    out = []
    for p in range(1, int((m + 1)*R2) + 1):
        for delta in (0.0, 0.5):
            f = frac(p/R2 + delta)
            if f < 1e-6: continue
            s = m + f
            qmax = int(math.floor(s*R2 - p + EPS))
            for q in range(max(p, qmax - 8), qmax + 1):
                if q < 1: continue
                t0 = (q/(2*R2) - s/2) % 0.5
                for t in {0.0, 0.25, -0.25, 0.5, round(t0, 12), round(t0 - 0.5, 12)}:
                    for order in (("TL", "BR", "BL", "TR"), ("BL", "TR", "TL", "BR")):
                        T = build(s, p, q, delta, t, order)
                        if T is None: continue
                        out.append((s, len(T), (p, q, delta, t, order)))
    return out

# ---------------------------------------------------------------------------
# Published best known sides for tilted packings with n <= 324, read in November 2026 from
# David Ellsworth's page (after Erich Friedman): https://kingbird.myphotos.cc/packing/squares_in_squares.html
# Any n <= 324 not listed has the plain grid as its best known packing.
PUBLISHED = """5 2.70711;10 3.70711;11 3.87708;17 4.67553;18 4.82288;19 4.88562;26 5.62132;27 5.70711;28 5.82444;29 5.93383;37 6.59862;38 6.70711;39 6.81072;40 6.82843;41 6.92669;50 7.57143;51 7.70080;52 7.70711;53 7.82288;54 7.84667;55 7.94577;65 8.53553;66 8.65685;67 8.70711;68 8.79880;69 8.82719;70 8.88167;71 8.94407;82 9.53553;83 9.63476;84 9.70711;85 9.74264;86 9.82288;87 9.83882;88 9.88815;89 9.94975;
101 10.53553;102 10.61139;103 10.70378;104 10.70711;105 10.80762;106 10.82298;107 10.84667;108 10.92592;109 10.94975;110 10.99679;
122 11.53553;123 11.60140;124 11.65685;125 11.70711;126 11.77474;127 11.82288;128 11.82509;129 11.88131;130 11.91119;131 11.95653;132 11.99137;
145 12.53553;146 12.60091;147 12.65685;148 12.65685;149 12.70711;150 12.77817;151 12.82288;152 12.83096;153 12.88167;154 12.93167;155 12.95845;156 12.98208;
170 13.53553;171 13.57143;172 13.61899;173 13.65685;174 13.70711;175 13.77817;176 13.82288;177 13.82303;178 13.84667;179 13.89534;180 13.93509;181 13.95691;182 13.97419;
197 14.53553;198 14.57143;199 14.61899;200 14.65685;201 14.70711;202 14.72792;203 14.77817;204 14.82288;205 14.82445;206 14.87222;207 14.89395;208 14.93777;209 14.95862;210 14.97413;
226 15.53553;227 15.57107;228 15.60902;229 15.65685;230 15.68293;231 15.70711;232 15.77817;233 15.77817;234 15.82288;235 15.82661;236 15.87608;237 15.91421;238 15.93966;239 15.95635;240 15.97556;241 15.99081;
257 16.53553;258 16.57107;259 16.60257;260 16.65685;261 16.68293;262 16.70711;263 16.74264;264 16.77817;265 16.77817;266 16.82306;267 16.84667;268 16.87931;269 16.90597;270 16.94062;271 16.95500;272 16.96972;273 16.98821;
290 17.53553;291 17.53553;292 17.60257;293 17.63415;294 17.65685;295 17.70711;296 17.70711;297 17.74106;298 17.77817;299 17.82288;300 17.82412;301 17.86889;302 17.88675;303 17.93126;304 17.94911;305 17.96066;306 17.96914;307 17.98272"""

def published():
    return {int(a): float(b) for a, b in (x.split() for x in PUBLISHED.replace("\n", "").split(";"))}

def close_table(tab, nmax):
    """Fewer tiles fit wherever more fit; and adding an L (one row and one column) grows the side by 1."""
    changed = True
    while changed:
        changed = False
        for n in range(1, nmax + 1):
            s = tab[n]; n2 = n + 2*int(math.floor(s + EPS)) + 1
            if n2 <= nmax and tab[n2] > s + 1 + EPS: tab[n2] = s + 1; changed = True
        lo = 1e9
        for n in range(nmax, 0, -1):
            if tab[n] < lo - 1e-12: lo = tab[n]
            elif tab[n] > lo + EPS: tab[n] = lo; changed = True
    return tab

def proven_floor(nmax):
    """sqrt(n), and for non-square n >= 8 the bound 1/2 + sqrt(n - floor(sqrt n) + 1/4) (Karakus 2026, as quoted by MathWorld)."""
    out = {}; f = 0.0
    for n in range(1, nmax + 1):
        r = math.isqrt(n); v = float(r) if r*r == n else math.sqrt(n)
        if r*r != n and n >= 8: v = max(v, 0.5 + math.sqrt(n - r + 0.25))
        f = max(f, v); out[n] = f
    return out

def main(lo=325, hi=1000, out="table_325_to_1000.csv"):
    import csv
    rec = published(); pub = lambda n: rec.get(n, float(math.ceil(math.sqrt(n))))
    lim = hi + 80
    # --- what the published records below 325 give when carried upward
    carry = close_table({n: (pub(n) if n <= 324 else float(math.ceil(math.sqrt(n)))) for n in range(1, hi + 1)}, hi)
    # --- one-row pattern: waste/m as a function of x/m, taken from blocks m = 15, 16, 17
    us = np.linspace(1/17, 1 + 1/17, 60); curves = []
    for m in (15, 16, 17):
        pts = [(x/m, (pub(m*m + x)**2 - (m*m + x))/m) for x in range(1, m + 2)]
        curves.append(np.interp(us, [a for a, b in pts], [b for a, b in pts]))
    cu = np.mean(curves, axis=0)
    def pattern(n):
        m = math.isqrt(n); x = n - m*m
        if x == 0: return float(m)
        if x > m + 1: return float(m + 1)
        return min(float(m + 1), math.sqrt(n + float(np.interp(x/m, us, cu))*m))
    # --- build band layouts for every block, keep the best count for each side
    top = {}
    for m in range(1, math.isqrt(hi) + 1):
        for s, c, rc in search_block(m):
            k = round(s, 9)
            if k not in top or c > top[k][1]: top[k] = (s, c, rc)
    best = {k*k: (float(k), ("grid", k), []) for k in range(1, math.isqrt(lim) + 2)}
    for s, c, rc in top.values():
        if c <= lim and (c not in best or s < best[c][0] - EPS): best[c] = (s, ("band", s, rc), [])
    changed = True
    while changed:
        changed = False; cur = None
        for n in range(lim, 0, -1):
            if n in best and (cur is None or best[n][0] < cur[1][0] - EPS): cur = (n, best[n])
            elif cur is not None and (n not in best or cur[1][0] < best[n][0] - EPS):
                best[n] = (cur[1][0], cur[1][1], cur[1][2] + [("drop", cur[0] - n)]); changed = True
        for n in sorted(best):
            s, b, ops = best[n]; n2 = n + 2*int(math.floor(s + EPS)) + 1
            if n2 <= lim and (n2 not in best or s + 1 < best[n2][0] - EPS):
                best[n2] = (s + 1, b, ops + [("L",)]); changed = True
    def realize(n):
        s, b, ops = best[n]
        if b[0] == "grid":
            k = b[1]; g = np.arange(k) + .5; X, Y = np.meshgrid(g, g); T = np.c_[X.ravel(), Y.ravel(), np.zeros(k*k)]; side = float(k)
        else:
            side = b[1]; p, q, d, t, order = b[2]; T = build(side, p, q, d, t, order)
        for op in ops:
            if op[0] == "L": T, side = add_L(T, side)
            elif op[1]: T = T[:len(T) - op[1]]
        return T, side
    def label(n):
        s, b, ops = best[n]
        nl = sum(1 for o in ops if o[0] == "L"); nd = sum(o[1] for o in ops if o[0] == "drop")
        if b[0] == "grid": base = f"plain grid {b[1]}x{b[1]}"
        else:
            p, q, d, t, order = b[2]; base = f"45deg band {p}x{q}" + (" shifted half" if d else "") + f" in side {b[1]:.4f}"
        return base + (f", +{nl} L" if nl else "") + (f", -{nd} tiles" if nd else "")
    # --- safety checks
    below = [n for n in range(1, 325) if best[n][0] < pub(n) - 2e-5]
    assert not below, ("built below a published record, so something is wrong", below)
    floor = proven_floor(hi); worst = 0.0
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["n", "as m^2+x", "proven floor", "predicted side", "prediction based on", "built and checked side", "how it was built", "built minus predicted", "built minus floor"])
        for n in range(lo, hi + 1):
            m = math.isqrt(n); x = n - m*m; grid = float(math.ceil(math.sqrt(n)))
            T, side = realize(n); assert len(T) == n
            worst = max(worst, check(T, side))
            pat = pattern(n); car = carry[n]
            if x == 0: pred, basis = float(m), "perfect square"
            elif car <= pat + EPS and car < grid - EPS: pred, basis = car, "published record below 325, carried up by adding L"
            elif pat < grid - EPS: pred, basis = pat, "one-row pattern fitted on blocks 15 to 17"
            else: pred, basis = grid, "plain grid expected"
            w.writerow([n, f"{m}^2+{x}", f"{floor[n]:.4f}", f"{pred:.4f}", basis, f"{side:.4f}", label(n), f"{side - pred:+.4f}", f"{side - floor[n]:.4f}"])
    assert worst < 1e-9, worst
    print("wrote", out, "- every layout checked, worst overlap", worst)

if __name__ == "__main__":
    main()
