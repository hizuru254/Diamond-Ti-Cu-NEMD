#!/usr/bin/env python3
"""Pool every flat-interface measurement to get the honest central value and spread."""
import numpy as np
vals = {
 "v001 seed1001": 7.292e-10, "v001 seed1002": 7.464e-10, "v001 seed1003": 7.229e-10,
 "v005 rebuilt":  6.931e-10, "v007 duplicate": 6.931e-10,
 "v008 flux3":    7.642e-10, "v009 flux2":     6.290e-10,
}
v = np.array(list(vals.values()))
print("individual (1e-10):", np.round(v*1e10, 3))
print(f"n={len(v)}  mean={v.mean():.4e}  std={v.std(ddof=1):.3e} "
      f"({100*v.std(ddof=1)/v.mean():.1f} %)  min={v.min():.3e}  max={v.max():.3e}")
print(f"sem={v.std(ddof=1)/np.sqrt(len(v)):.3e}")
