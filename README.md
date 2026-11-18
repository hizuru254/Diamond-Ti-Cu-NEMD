# Si/Ti/Cu/Ti/Si 界面热阻（NEP89 势）

**目标**：用 NEMD 算 Si/Ti/Cu/Ti/Si 叠层的界面热阻（TBR），并做 Au/Ti 的替代研究——
即把原项目的 Au 中间层换成 Cu，用真实多元素势（NEP89）替代 Tersoff+EAM+LJ 混合。

**来源**：本仓库继承自 `hizuru254/Si-Ti-Au-NEMD` 的流程与工具链（协议 v1.6），
只替换材料与势函数。

---

## 0. 为什么换材料 / 换势

原项目（Si/Ti/Au/Ti/Si）的三个硬伤，都是势函数造成的：

1. Ti–Si、Ti–Au 只有手挑的 LJ 交叉项（ε=0.2–0.3 eV、σ=2.0 Å）。强 → 50 ps 内
   互混放热 2206 eV；弱 → 层间不粘、Si 板膨胀 15 Å。**没有物理合理的中间值**。
2. 界面持续弛豫 → TBR 单调上涨（160 ps 内 +73%），`time_converged ≤5%` 不可达。
3. 提取方法之间差 2.6–5 倍（窗口外推偏低），界面宽 ~10 Å 无法定义。

现在改用 **NEP89**（89 元素通用神经演化势，含 Si / Ti / Cu）：界面化学由第一性
原理数据训练，不再需要用 LJ 硬凑。

## 1. 势函数

| 项 | 值 |
|---|---|
| 模型 | `potentials/nep89.txt`（14 MB，976,331 参数） |
| 来源 | GPUMD 官方 `potentials/nep/nep89_20250409/` |
| 获取 | `bash tools/get_nep89.sh` |
| 引用 | Liang et al., *NEP89: Universal neuroevolution potential...*, Nature Comput. Sci. (2026) |
| LAMMPS | 需要 `USER-NEP` 包（本机已编译：`/home/cc/lammps/src/lmp_mpi`） |

编译方法（25 秒）：

```bash
git clone --depth 1 https://github.com/brucefan1983/NEP_CPU.git /tmp/NEP_CPU
cd /home/cc/lammps/src
cp -r /tmp/NEP_CPU/interface/lammps/USER-NEP .
make yes-USER-NEP && make -j32 mpi
```

## 2. 已实测的基线（不是猜的）

| 项 | 结果 |
|---|---|
| NEP89 元素覆盖 | 含 Si、Ti、Cu、Au ✓ |
| LAMMPS+NEP 编译 | 25 s，`lmp_mpi` 可跑 ✓ |
| Cu 平衡晶格 | **a₀ = 3.60 Å**（实验 3.615，差 0.4%） |
| 性能 | 8000 原子 / 8 ranks = 0.0935 s/步 → **3 万原子 / 24 ranks ≈ 8.5 步/秒** |
| 150 ps 单点成本 | ≈ **5 小时**（24 ranks） |

## 3. 可信度要求（必须做完才报数）

1. **体相标定**：Si、Ti、Cu 的 a₀ / c / 内聚能 / 弹性常数 vs 实验或 DFT。
   写进 `potentials.lock.json` 的 `bulk_validation`（现在是 WAIVED，不可接受）。
2. **单界面优先**：先做 Ti/Cu 与 Si/Ti 两层算例（`tools/build_pair.py`），
   两侧各 ≥60 Å，ΔT 用两个平台之差，不依赖外推。
3. **串联验证**：3 层 `R(Si/Ti)+R(Ti/Cu)` 与直接测的 3 层值对照。
4. **5 层对照**：`2×(Si/Ti + Ti/Cu)` 与 5 层直测在误差内一致。
5. **文献对照**：TBR 量级与已发表值比较（Si/Ti、Cu/Ti 各有 TDTR/MD 数据）。

## 4. 状态

- [x] 仓库骨架、工具链、NEP89 接入
- [ ] **体相标定**（下一步，`tools/calibrate.py`）
- [ ] 单界面 Ti/Cu
- [ ] 单界面 Si/Ti
- [ ] 3 层串联验证
- [ ] 5 层锚点
- [ ] 图案扫描（槽宽/深度/周期，Au→Cu 的界面）

## 5. 继承自原项目、已经验证过的流程

- 协议 v1.6：run 内退火 300→700→300 K → 平衡 → 通量热身 80 ps（不计量）→ 计量 150 ps
- 三种 TBR 估计量：`tbr_plane.py`（首选）/ `tbr_global_fit.py` / `tbr_uncertainty.py`
- 并行调度 `scan_parallel.py`（**并发 rank 总数 ≤48**，超过会慢 40 倍且无报错）
- 不可覆盖的 run 目录 + manifest 哈希 + ledger

详见 `docs/PLAN.md` 与继承的 `docs/PROTOCOL_V16.md`。
