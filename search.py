"""Shake-and-squeeze search: try to fit n unit squares into a smaller square.

Start from a known layout, disturb it, squeeze it again, and keep the result if the box got smaller.
Every result that is saved has been checked tile against tile (see band_table.check).

    python search.py --n 419 --minutes 300 --seed 1 --out results

This is a slow, simple search written for clarity. It is not expected to match the long-running
programs behind the published records. It is a first tool to see where the box can shrink.
"""
import argparse, csv, math, os, time
# one thread per search: several searches run side by side, and extra threads only get in each other's way
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")
import numpy as np
from scipy.optimize import minimize
import band_table as bt


def pairs(x, y, r=1.9):
    d2 = (x[:, None] - x[None, :])**2 + (y[:, None] - y[None, :])**2
    return np.nonzero(np.triu(d2 < r*r, 1))


def energy(v, n, s, I, J):
    """Sum of squared overlaps (tile-tile and tile-wall), its gradient, and its slope in the side s."""
    x = v[:n]; y = v[n:2*n]; t = v[2*n:3*n]
    c = np.cos(t); sn = np.sin(t)
    g = np.zeros(3*n); gs = 0.0; E = 0.0
    if len(I):
        dx = x[J] - x[I]; dy = y[J] - y[I]
        ci = c[I]; si = sn[I]; cj = c[J]; sj = sn[J]
        cphi = ci*cj + si*sj; sphi = si*cj - ci*sj
        hh = .5 + .5*(np.abs(cphi) + np.abs(sphi))
        P = np.stack([dx*ci + dy*si, -dx*si + dy*ci, dx*cj + dy*sj, -dx*sj + dy*cj])
        O = hh - np.abs(P)
        k = np.argmin(O, axis=0); d = O[k, np.arange(len(I))]
        act = d > 0
        if act.any():
            Ia = I[act]; Ja = J[act]; da = d[act]; ka = k[act]; Pa = P[:, act]; idx = np.arange(len(da))
            sg = np.sign(Pa[ka, idx])
            cia, sia, cja, sja = ci[act], si[act], cj[act], sj[act]
            ax = np.choose(ka, [cia, -sia, cja, -sja]); ay = np.choose(ka, [sia, cia, sja, cja])
            E += float(np.sum(da*da)); w = 2*da
            fx = w*sg*ax; fy = w*sg*ay
            g[:n] += np.bincount(Ia, fx, n) - np.bincount(Ja, fx, n)
            g[n:2*n] += np.bincount(Ia, fy, n) - np.bincount(Ja, fy, n)
            dhh = .5*(-np.sign(cphi[act])*sphi[act] + np.sign(sphi[act])*cphi[act])
            z = np.zeros_like(da)
            dpi = np.choose(ka, [Pa[1], -Pa[0], z, z]); dpj = np.choose(ka, [z, z, Pa[3], -Pa[2]])
            g[2*n:] += np.bincount(Ia, w*(dhh - sg*dpi), n) + np.bincount(Ja, w*(-dhh - sg*dpj), n)
    h = .5*(np.abs(c) + np.abs(sn)); dh = .5*(-np.sign(c)*sn + np.sign(sn)*c)
    for pos, off in ((x, 0), (y, n)):
        w1 = np.maximum(h - pos, 0); w2 = np.maximum(pos + h - s, 0)
        E += float(np.sum(w1*w1) + np.sum(w2*w2))
        g[off:off+n] += -2*w1 + 2*w2
        g[2*n:] += 2*(w1 + w2)*dh
        gs += -2*float(np.sum(w2))
    return E, g, gs


def worst_overlap(v, n, s):
    T = np.c_[v[:n], v[n:2*n], np.degrees(v[2*n:3*n])]
    return bt.check(T, s)


def relax(v, n, s, iters=1500):
    I, J = pairs(v[:n], v[n:2*n])
    f = lambda q: energy(q, n, s, I, J)[:2]
    return minimize(f, v, jac=True, method="L-BFGS-B", options=dict(maxiter=iters, ftol=1e-16, gtol=1e-12, maxcor=20)).x


def joint(v, n, s, mu, iters=2500):
    I, J = pairs(v[:n], v[n:2*n])
    def f(q):
        E, g, gs = energy(q[:-1], n, q[-1], I, J)
        return q[-1] + mu*E, np.append(mu*g, 1 + mu*gs)
    r = minimize(f, np.append(v, s), jac=True, method="L-BFGS-B", options=dict(maxiter=iters, ftol=1e-16, gtol=1e-12, maxcor=20)).x
    return r[:-1], float(r[-1])


def squeeze(v, n, s, stages=(1e2, 1e4, 1e6, 1e8), reps=2):
    """Shrink the box until the tiles jam. Returns a strictly valid layout, its side, and the overlap left (0)."""
    v = relax(v, n, s)
    # pull the box in to touch the tiles, then a little further, so the squeeze has something to push on
    c = np.cos(v[2*n:]); sn = np.sin(v[2*n:]); h = .5*(np.abs(c) + np.abs(sn))
    v = v.copy(); v[:n] -= (v[:n] - h).min(); v[n:2*n] -= (v[n:2*n] - h).min()
    s = max((v[:n] + h).max(), (v[n:2*n] + h).max()) - 1e-4
    for mu in stages:
        for _ in range(reps):
            v, s = joint(v, n, s, mu)
    # the penalty method leaves tiny overlaps: grow everything by a hair until none is left
    for _ in range(8):
        o = worst_overlap(v, n, s)
        if o <= 1e-12: break
        lam = 1 + 3*o + 1e-12
        v = v.copy(); v[:2*n] *= lam; s *= lam
    return v, s, worst_overlap(v, n, s)


def shake(v, n, s, rng):
    """One random disturbance. Returns the new layout and how much extra room to give it."""
    v = v.copy(); x = v[:n]; y = v[n:2*n]; t = v[2*n:]
    kind = rng.choice(["jiggle", "turn", "slide", "hop", "big_turn"], p=[.25, .3, .15, .2, .1])
    if kind == "jiggle":
        a = rng.choice([0.01, 0.03, 0.08])
        x += rng.normal(0, a, n); y += rng.normal(0, a, n); t += rng.normal(0, a, n)
        return v, 3*a + 0.02
    i = rng.integers(n); R = rng.uniform(1.2, 4.5) if kind != "big_turn" else rng.uniform(3, 8)
    grp = np.hypot(x - x[i], y - y[i]) < R
    if kind in ("turn", "big_turn"):
        r = math.radians(rng.normal(0, 6 if kind == "turn" else 12)); cx, cy = x[grp].mean(), y[grp].mean()
        dx, dy = x[grp] - cx, y[grp] - cy
        x[grp] = cx + dx*math.cos(r) - dy*math.sin(r); y[grp] = cy + dx*math.sin(r) + dy*math.cos(r); t[grp] += r
        return v, 0.08 + abs(r)*R*0.5
    if kind == "slide":
        a = rng.uniform(0, 2*math.pi); d = rng.uniform(0.05, 0.4)
        x[grp] += d*math.cos(a); y[grp] += d*math.sin(a)
        return v, 0.05 + d
    # hop: move one tile into the emptiest spot found
    pts = rng.uniform(0.5, s - 0.5, (400, 2))
    dist = np.min(np.hypot(pts[:, None, 0] - x[None, :], pts[:, None, 1] - y[None, :]), axis=1)
    j = int(np.argmax(dist)); x[i], y[i] = pts[j]; t[i] = rng.uniform(-math.pi/4, math.pi/4)
    return v, 0.15


def save(path, n, s, v, note):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        fh.write(f"# n {n}\n# side {s:.9f}\n# {note}\n# tile,x,y,angle_deg\n")
        w = csv.writer(fh)
        for k in range(n):
            w.writerow([k + 1, f"{v[k]:.9f}", f"{v[n+k]:.9f}", f"{math.degrees(v[2*n+k]):.7f}"])


def load(path):
    side = None; rows = []
    for line in open(path):
        if line.startswith("# side"): side = float(line.split()[2])
        elif not line.startswith("#") and line.strip():
            p = line.split(","); rows.append((float(p[1]), float(p[2]), math.radians(float(p[3]))))
    a = np.array(rows); return np.concatenate([a[:, 0], a[:, 1], a[:, 2]]), side, len(rows)


def offcentre_start(n):
    """Smallest 45-degree band layout that sits half a tile off centre and holds at least n tiles."""
    m = math.isqrt(n); best = None
    for s, c, rc in bt.search_block(m):
        if rc[2] == 0.5 and c >= n and (best is None or s < best[0] - 1e-9 or (abs(s - best[0]) < 1e-9 and c < best[1])):
            best = (s, c, rc)
    if best is None: return None
    s, c, rc = best; T = bt.build(s, *rc)[:n]
    return T, s, f"off-centre 45deg band {rc[0]}x{rc[1]} in side {s:.4f}, -{c - n} tiles"


def run(n, minutes, seed, out, start=None, log=print, offcentre=False):
    rng = np.random.default_rng(seed); t_end = time.time() + 60*minutes
    oc = offcentre_start(n) if offcentre else None
    if start and os.path.exists(start):
        v, s, k = load(start); assert k == n, "start file has the wrong number of tiles"; origin = f"start file {start}"
    elif oc is not None:
        T, s, origin = oc; assert bt.check(T, s) < 1e-9
        v = np.concatenate([T[:, 0], T[:, 1], np.radians(T[:, 2])])
    else:
        best_tab = bt.build_best(max(n, 30))
        T, s = bt.realize(best_tab, n); origin = bt.label(best_tab, n)
        v = np.concatenate([T[:, 0], T[:, 1], np.radians(T[:, 2])])
    s0 = s
    cur_v, cur_s, o = squeeze(v, n, s + 1e-6)
    if o > 1e-11 or cur_s > s0: cur_v, cur_s = v, s0      # never start worse than the checked start
    best_v, best_s = cur_v, cur_s
    path = os.path.join(out, f"{n}.csv")
    save(path, n, best_s, best_v, f"seed {seed}; started from: {origin} (side {s0:.6f})")
    log(f"n={n} start {s0:.6f} ({origin}); after first squeeze {best_s:.6f}")
    steps = 0; temp = 0.004
    while time.time() < t_end:
        steps += 1
        cand, room = shake(cur_v, n, cur_s, rng)
        try:
            # quick, rough squeeze first; only a promising result gets the slow, exact squeeze
            cv, cs, o = squeeze(cand, n, cur_s + room, stages=(1e3, 1e5), reps=1)
            if o <= 1e-11 and cs < cur_s + 0.004:
                cv, cs, o = squeeze(cv, n, cs, stages=(1e6, 1e8))
        except Exception as e:                      # a failed step is skipped, never saved
            log("step failed:", e); continue
        if o > 1e-11: continue
        if cs < cur_s - 1e-9 or rng.random() < math.exp(-(cs - cur_s)/temp):
            cur_v, cur_s = cv, cs
        if cs < best_s - 1e-8:
            best_v, best_s = cv, cs
            save(path, n, best_s, best_v, f"seed {seed}; started from: {origin} (side {s0:.6f}); step {steps}")
            log(f"n={n} step {steps}: new best {best_s:.6f}")
        temp = max(0.0005, temp*0.999)
    log(f"n={n} done: {steps} steps, start {s0:.6f}, best {best_s:.6f}")
    return best_s


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--minutes", type=float, default=10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="results")
    ap.add_argument("--start", default=None, help="optional layout file to start from (same format as the results)")
    ap.add_argument("--offcentre", action="store_true", help="start from a band that sits half a tile off centre")
    a = ap.parse_args()
    run(a.n, a.minutes, a.seed, a.out, a.start, offcentre=a.offcentre)
