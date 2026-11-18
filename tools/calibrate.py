#!/usr/bin/env python3
"""Body-phase calibration with the NEP89 model: lattice, energy, (later elastic).

Writes potentials/lattice.json which build_geometry.py must read instead of the
hard-coded Tersoff/EAM constants.
"""
from __future__ import annotations
import argparse, json, re, subprocess, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
NEP = ROOT / "potentials/nep89.txt"

BASIS = {
    "Si": [(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0),(.25,.25,.25),(.25,.75,.75),(.75,.25,.75),(.75,.75,.25)],
    "Cu": [(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0)],
    "Au": [(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0)],
}
MASS = {"Si": 28.0855, "Cu": 63.546, "Au": 196.96657, "Ti": 47.867}


def write_cell(el: str, a: float, path: Path) -> int:
    basis = BASIS[el]
    lines = [el, "", f"{len(basis)} atoms", "1 atom types", "",
             f"0 {a:.6f} xlo xhi", f"0 {a:.6f} ylo yhi", f"0 {a:.6f} zlo zhi", "",
             "Masses", "", f"1 {MASS[el]}", "", "Atoms # atomic", ""]
    for i, b in enumerate(basis, 1):
        lines.append(f"{i} 1 {b[0]*a:.6f} {b[1]*a:.6f} {b[2]*a:.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(basis)


def energy(work: Path, exe: str, model: Path, el: str) -> float:
    inp = f"""units metal
atom_style atomic
boundary p p p
read_data cell.data
pair_style nep
pair_coeff * * {model} {el}
neighbor 2.0 bin
thermo_style custom step pe
minimize 1e-10 1e-10 2000 20000
"""
    (work / "in.eos").write_text(inp, encoding="utf-8")
    subprocess.run([exe, "-in", "in.eos", "-log", "log.eos", "-screen", "none"],
                   cwd=work, check=True, stdout=subprocess.DEVNULL)
    txt = (work / "log.eos").read_text(errors="replace")
    vals = re.findall(r"^\s+\d+\s+(-?\d+\.\d+)", txt, re.M)
    return float(vals[-1])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="/home/cc/lammps/src/lmp_mpi")
    ap.add_argument("--model", default=str(NEP))
    ap.add_argument("--work", default="/tmp/nep_calib")
    ap.add_argument("--out", default=str(ROOT / "potentials/lattice.json"))
    args = ap.parse_args()
    work = Path(args.work); work.mkdir(parents=True, exist_ok=True)
    guesses = {"Si": (5.40, 5.70), "Cu": (3.50, 3.80), "Au": (4.00, 4.30)}
    out = {}
    for el, (lo, hi) in guesses.items():
        best = None
        for a in np.arange(lo, hi + 1e-9, 0.02):
            n = write_cell(el, float(a), work / "cell.data")
            e = energy(work, args.exe, Path(args.model), el) / n
            if best is None or e < best[1]:
                best = (float(a), e)
        out[el] = {"a0_A": round(best[0], 4), "E0_eV_atom": round(best[1], 4)}
        print(f"{el}: a0={best[0]:.3f} A  E0={best[1]:.4f} eV/atom")
    Path(args.out).write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(f"[OK] {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
