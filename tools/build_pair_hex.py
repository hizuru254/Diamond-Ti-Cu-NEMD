#!/usr/bin/env python3
"""Cu(111)/Ti(0001) single-interface cell, exact 7:8 commensuration, numpy only."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

A_TI, C_TI = 2.950, 4.602
A0_CU = 3.650
A_CU = A0_CU/np.sqrt(2.0)     # 2.5810  NN on Cu(111)
D_CU = A0_CU/np.sqrt(3.0)     # 2.1075  interlayer spacing

def layer(nx, ny, a, shift):
    pts = []
    for i in range(nx):
        for j in range(ny):
            for bx, by in ((0.0,0.0), (0.5,0.5)):
                pts.append((((i+bx)*a + shift[0]) % (nx*a),
                            ((j+by)*a*np.sqrt(3) + shift[1]) % (ny*a*np.sqrt(3))))
    return np.array(pts)

def slab(kind, nz, a, dz, z0):
    nx = ny = 7 if kind == "ti" else 8
    out = []
    for k in range(nz):
        if kind == "ti":
            sh = (a/2.0, a/(2.0*np.sqrt(3))) if k % 2 else (0.0, 0.0)
        else:
            s = k % 3
            sh = (s*a/2.0, s*a/(2.0*np.sqrt(3)))
        p = layer(nx, ny, a, sh)
        out.append(np.hstack([p, np.full((len(p),1), z0 + k*dz)]))
    return np.vstack(out)

def min_dist(x, y, Lx, Ly, chunk=512, same=False):
    """Minimum periodic distance between two point sets (self-pairs excluded)."""
    m = np.inf
    for i in range(0, len(x), chunk):
        dx = np.abs(x[i:i+chunk,None,0] - y[None,:,0]); dx = np.minimum(dx, Lx-dx)
        dy = np.abs(x[i:i+chunk,None,1] - y[None,:,1]); dy = np.minimum(dy, Ly-dy)
        dz = x[i:i+chunk,None,2] - y[None,:,2]
        d = np.sqrt(dx*dx + dy*dy + dz*dz)
        if same:
            for k in range(d.shape[0]):
                d[k, i+k] = np.inf
        m = min(m, d.min())
    return m


def _groups(ti, cu, hot_thick=0.0):
    """4 Ti layers (4*98=392 atoms) against 3 Cu layers (3*128=384): 2 % apart."""
    zt, zc = ti[:,2], cu[:,2]
    h_lo = float(zt.min()) + 8.0
    h_hi = h_lo + 4*(C_TI/2) + 0.6
    return {"hot": {"z_min": h_lo, "z_max": h_hi},
            "cold": {"z_min": float(zc.max()) - 3*D_CU - 0.6, "z_max": float(zc.max())}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--thick-ti", type=float, default=60.0)
    ap.add_argument("--thick-cu", type=float, default=60.0)
    ap.add_argument("--vacuum", type=float, default=15.0)
    ap.add_argument("--gap", type=float, default=2.45, help="target Ti-Cu NN across the interface")
    ap.add_argument("--bin", type=float, default=2.0)
    ap.add_argument("--groove-w", type=float, default=0.0, help="groove width along x (A)")
    ap.add_argument("--groove-d", type=float, default=0.0, help="groove depth into Cu (A)")
    ap.add_argument("--groove-p", type=float, default=0.0, help="groove period along x (A)")
    ap.add_argument("--alloy", type=float, default=0.0,
                    help="thickness (A) of a graded Ti-Cu alloy layer at the interface")
    ap.add_argument("--alloy-seed", type=int, default=7)
    ap.add_argument("--interlock", type=float, default=0.0,
                    help="interlocked corrugation height H (A); requires --groove-w and --groove-p")
    ap.add_argument("--out-structure", required=True)
    ap.add_argument("--out-geom", required=True)
    a = ap.parse_args()

    lx_ti, ly_ti = 7*A_TI, 7*A_TI*np.sqrt(3)
    lx_cu, ly_cu = 8*A_CU, 8*A_CU*np.sqrt(3)
    print(f"[commensuration] Ti {lx_ti:.4f}x{ly_ti:.4f}  Cu {lx_cu:.4f}x{ly_cu:.4f} "
          f" mismatch {100*(lx_cu-lx_ti)/lx_ti:+.3f}%")
    if abs(lx_cu-lx_ti) > 0.02*lx_ti or abs(ly_cu-ly_ti) > 0.02*ly_ti:
        print("[ERROR] in-plane mismatch too large"); return 1
    L, Ly = lx_ti, ly_ti

    nz_ti = int(round(a.thick_ti/(C_TI/2)))
    nz_cu = int(round(a.thick_cu/D_CU))
    ti = slab("ti", nz_ti, A_TI, C_TI/2, 0.0)
    ti[:,2] -= ti[:,2].max()                       # Ti top layer at z = 0
    z_top_ti = 0.0

    # --- optional interlocked corrugation ----------------------------------
    H = a.interlock
    if H > 0 and a.groove_w > 0 and a.groove_p > 0:
        phase = np.mod(ti[:,0], a.groove_p) - 0.5*a.groove_p
        ridge = np.abs(phase) < 0.5*a.groove_w
        # cut trenches into Ti outside the ridge columns
        cut = (~ridge) & (ti[:,2] > z_top_ti - H)
        print(f"[interlock] H={H} W={a.groove_w} P={a.groove_p}: cutting {int(cut.sum())} Ti atoms")
        ti = ti[~cut]
        z_cu_start = z_top_ti - H + a.gap
    else:
        ridge = None
        z_cu_start = z_top_ti + a.gap

    cu = slab("cu", nz_cu + int(np.ceil((z_cu_start)/D_CU)) + 2, A_CU, D_CU, z_cu_start)
    cu = cu[(cu[:,2] >= z_cu_start) & (cu[:,2] < z_cu_start + a.thick_cu + 3*D_CU)]
    if ridge is not None:
        phase_c = np.mod(cu[:,0], a.groove_p) - 0.5*a.groove_p
        ridge_c = np.abs(phase_c) < 0.5*a.groove_w
        clash = ridge_c & (cu[:,2] < z_top_ti + 0.5*a.gap)
        print(f"[interlock] removing {int(clash.sum())} Cu atoms that would overlap Ti ridges")
        cu = cu[~clash]
    # --- optional periodic grooves cut into the Cu side --------------------
    if a.groove_w > 0 and a.groove_d > 0 and a.groove_p > 0:
        z_int = 0.5*(float(ti[:,2].max()) + float(cu[:,2].min()))
        phase = np.mod(cu[:,0], a.groove_p) - 0.5*a.groove_p
        in_groove = (np.abs(phase) < 0.5*a.groove_w) & (cu[:,2] < z_int + a.groove_d)
        print(f"[groove] W={a.groove_w} D={a.groove_d} P={a.groove_p}: "
              f"removing {int(in_groove.sum())} Cu atoms")
        cu = cu[~in_groove]

    # drop Cu atoms that sit closer than 1.8 A to any Ti atom (interface overlap)
    keep = np.ones(len(cu), bool)
    for i in range(0, len(cu), 1024):
        dx = np.abs(cu[i:i+1024,None,0] - ti[None,:,0]); dx = np.minimum(dx, L-dx)
        dy = np.abs(cu[i:i+1024,None,1] - ti[None,:,1]); dy = np.minimum(dy, Ly-dy)
        dz = np.abs(cu[i:i+1024,None,2] - ti[None,:,2])
        close = ((dx**2 + dy**2 + dz**2) < 2.3**2).any(axis=1)
        keep[i:i+1024] = ~close
    n_drop = int((~keep).sum())
    cu = cu[keep]

    atoms = np.vstack([ti, cu])
    types = np.concatenate([np.full(len(ti),2), np.full(len(cu),3)])

    # --- optional graded Ti-Cu alloy layer ---------------------------------
    if a.alloy > 0:
        rng = np.random.default_rng(a.alloy_seed)
        z_all = atoms[:,2]
        zi_mid = 0.5*(z_all[types == 2].max() + z_all[types == 3].min())
        lo, hi = zi_mid - 0.5*a.alloy, zi_mid + 0.5*a.alloy
        inreg = (z_all >= lo) & (z_all <= hi)
        frac = np.clip((z_all - lo)/max(a.alloy, 1e-9), 0.0, 1.0)   # 0 = Ti side
        u = rng.random(len(atoms))
        # Ti atom becomes Cu with probability frac; Cu atom becomes Ti with 1-frac
        to_cu = inreg & (types == 2) & (u < frac)
        to_ti = inreg & (types == 3) & (u > frac)
        n_cu = int(to_cu.sum()); n_ti = int(to_ti.sum())
        types = types.copy()
        types[to_cu] = 3
        types[to_ti] = 2
        print(f"[alloy] graded layer {a.alloy} A: {n_cu} Ti->Cu, {n_ti} Cu->Ti "
              f"({int(inreg.sum())} atoms in the region)")
    zmin = atoms[:,2].min() - a.vacuum
    zmax = atoms[:,2].max() + a.vacuum

    # --- self checks -------------------------------------------------------
    d_ii = min_dist(ti, ti, L, Ly, same=True)
    d_cc = min_dist(cu, cu, L, Ly, same=True)
    d_ic = min_dist(ti, cu, L, Ly)
    print(f"[geometry] Ti {len(ti)}  Cu {len(cu)} (dropped {n_drop})  total {len(atoms)}")
    print(f"[check] min Ti-Ti {d_ii:.3f} A | min Cu-Cu {d_cc:.3f} A | min Ti-Cu {d_ic:.3f} A")
    ok = d_ii > 1.6 and d_cc > 1.6 and d_ic > 1.6
    print(f"[check] {'PASS' if ok else 'FAIL'}: no pair below 1.6 A")
    if not ok: return 1

    # --- data file ---------------------------------------------------------
    lines = ["Cu/Ti interface cell", "", f"{len(atoms)} atoms", "3 atom types", "",
             f"0 {L:.6f} xlo xhi", f"0 {Ly:.6f} ylo yhi",
             f"{zmin:.6f} {zmax:.6f} zlo zhi", "",
             "Masses", "", "1 12.011", "2 47.867", "3 63.546", "",
             "Atoms # atomic", ""]
    for i,(p,t) in enumerate(zip(atoms, types), 1):
        lines.append(f"{i} {t} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
    Path(a.out_structure).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_structure).write_text("\n".join(lines)+"\n", encoding="utf-8")

    zi = 0.5*(float(ti[:,2].max()) + float(cu[:,2].min()))
    geom = {
      "schema_version":"1.0","status":"APPROVED","geometry_id":"PAIR_CuTi_hex",
      "note":"Cu(111)/Ti(0001) single interface, 7:8 commensuration, no applied strain",
      "box":{"lx":L,"ly":Ly,"z_min":zmin,"z_max":zmax,"boundary":"p p p"},
      "layers":[{"id":"TI_L","material":"Ti","z_min":float(ti[:,2].min()),"z_max":zi},
                {"id":"CU_R","material":"Cu","z_min":zi,"z_max":float(cu[:,2].max())}],
      "interfaces":[{"id":"TI_CU","z":zi,
                     "left":{"z_min":zi-16.0,"z_max":zi-2.0,"material":"Ti"},
                     "right":{"z_min":zi+2.0,"z_max":zi+16.0,"material":"Cu"}}],
      "groups":_groups(ti, cu, 12.0),
      "chunk_bin_A":a.bin,"expected_atoms":{"total":int(len(atoms))}}
    Path(a.out_geom).write_text(json.dumps(geom, indent=2, ensure_ascii=False)+"\n",
                                encoding="utf-8", newline="\n")
    print(f"[geom] box {L:.2f} x {Ly:.2f} x {zmax-zmin:.2f}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
