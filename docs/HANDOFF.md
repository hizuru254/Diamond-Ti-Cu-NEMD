# HANDOFF — 当前状态与下一步（2026-11-18）

> 换新对话先读这个。它记录**已验证的事实**、**我自己犯过的错**、以及下一步的确切动作。

## 1. 材料与势

- 材料：**C(金刚石)/Ti/Cu/Ti/C**（原 Si/Ti/Au 换成 金刚石/Ti/Cu）
- 势：`potentials/CuTiC_NEP.txt`（来自 wugai-whu/NEP-Cu-X-diamond-thermal-conductivity，
  元素 Cu C Cr Mo Ti W，24k 参数；由旧 `nep4` 格式转成 `format_version 1`）
- 运行环境：`/home/cc/lammps/src/lmp_mpi`（源码编译，带 `USER-NEP`，25 秒编完）

## 2. 已验证（用立方构造，无歧义）

| 元素 | a₀ (NEP) | 实验 | 偏差 | E₀ (eV/atom) | 实验 |
|---|---|---|---|---|---|
| Cu (fcc) | **3.650 Å** | 3.615 | +1.0% | −3.711 | −3.49 |
| C (diamond) | **3.550 Å** | 3.567 | −0.5% | −9.084 | −7.37（参考偏移） |

两者都是**平滑凸极小**（完整 E(a) 表见 git 历史里 calibrate_cutic.py 的输出）。

## 3. 未完成：hcp Ti 的标定

脚本 `tools/calib_ti_hex.py` 用**精确三斜六方胞**跑网格：
a₁=(a,0,0)、a₂=(−a/2,√3a/2,0)、c，原子 (0,0,0)+(1/3,2/3,1/2)，LAMMPS 里 `xy=−a/2`。
判据：极小点应接近 a=2.951 Å、c/a=1.587，且能量面**凸**、无孤立跳变。

## 4. 我犯过的两个错（后人别重踩）

1. `fix box/relax aniso 1e9 vmax …` —— `1e9` 是**目标压力（bar）**，等于 10⁹ bar 压缩，
   盒子必塌、邻居表溢出。正确写法：`fix box/relax aniso 0.0 vmax 0.0002`。
2. hcp 正交胞的 basis 写错（用了 (2/3,1/3,1/2) 或 (0,1/3,1/2)）。**正交胞 (a×a√3×c)
   不能表示 hcp**；hcp 必须用六方/三斜胞，或 2 原子正交胞。用错 basis 会得到
   "能量面乱跳、极小跑到 c/a=1.414" 的假象。
3. 教训：标定必须**打印完整能量面**并检查凸性，不能只取 min。

## 5. 下一步（按顺序，不要跳）

1. 完成 Ti 标定（`tools/calib_ti_hex.py`）→ 把 a、c 写进 `potentials/lattice.json`。
2. 改 `tools/build_geometry.py` 为 C/Ti/Cu/Ti/C：面内取
   L = n·a_C（金刚石）且尽量与 Cu 共格；每层厚度 ≥40 Å；热源/冷源内移 10 Å。
3. `preflight` 到 READY（注意 hot/cold 原子数差异 ≤12% 是 warning，不是 error）。
4. **单界面 Ti/Cu 50 ps 稳定性测试**（`tools/build_pair.py` 要先支持 C/Ti）：
   看 PE 是否走平、Ti 是否互混。这一步不过，不要往下走。
5. 单界面 C/Ti 同理。
6. 3 层 C/Ti/Cu 验证串联可加性；7. 5 层锚点；8. 图案扫描。

## 6. 关键数字（省得重测）

- LAMMPS+USER-NEP 编译：25 s（`make yes-USER-NEP && make -j32 mpi`）
- NEP89（14 MB，976k 参数）性能：8000 原子 / 8 ranks = 0.0935 s/步
- 小 NEP（24k 参数）应快 10 倍以上，需实测
- 并发 rank 上限 **48**（超过慢 40 倍且无报错）
- 协议 v1.6：退火 300→700→300 K（65 ps）→ 平衡 10 ps → 通量热身 80 ps（不计量）→ 计量 150 ps
- 三种 TBR 估计量：`tbr_plane.py`（首选，无界面假设）/ `tbr_global_fit.py` / `tbr_uncertainty.py`
  （窗口外推 `postprocess.py` 在本类体系偏低 2.6–5 倍，只能作参考）

## 7. 可信度闸门（不满足就不报数）

1. Ti/Si 体相：a、c/a、E₀、弹性常数 vs 实验
2. 单界面 PE 走平（dPE/dt < 0.05 eV/ps）
3. 时间收敛 ≤5%（或给出统计+系统误差）
4. 串联可加性验证（3 层 vs 单界面之和）
5. 与文献 TBR 量级对照

## 8. 性能坑（2026-11-18 实测）

- `/home/cc/lammps/src/lmp_mpi` 是 **LAMMPS 2 Apr 2025 开发版**（不是 29Sep2021；
  那是 `/usr/bin/lmp` 发行版）。`config/params.json` 的 `version_marker` 要改成
  `"LAMMPS (2 Apr 2025"`，否则 preflight 的版本检查会 FAIL。
- **NEP 的固定开销很大**：54 原子单次力评估约 1 s（NEP89，976k 参数），
  所以 `minimize` 默认几千次迭代 → 单个点要几分钟。标定脚本必须：
  `minimize 1e-8 1e-8 200 2000`（少迭代）或干脆只算刚性胞能量（不 minimize）。
- 小 NEP（24k 参数）应快 ~40 倍，但仍需实测；标定时先跑单点计时再开网格。
