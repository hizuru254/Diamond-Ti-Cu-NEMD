# Diamond(C)/Ti/Cu/Ti/C 界面热阻：NEMD 图案化研究

用非平衡分子动力学（NEMD）计算 **金刚石(C)/Ti/Cu/Ti/C** 叠层的界面热阻（TBR），
并用同一套可复现流程寻找**界面图案化的最优几何**，给出可直接交付实验参考的
热阻-几何数据。

[![status](https://img.shields.io/badge/gates-G0%20READY-brightgreen)]()
[![potential](https://img.shields.io/badge/potential-NEP%20Cu%2FTi%2FC-blue)]()
[![baseline](https://img.shields.io/badge/TBR-7.33e--10%20m2K%2FW-orange)]()

---

## 1. 为什么是这个体系、这个势

原项目做 Si/Ti/Au/Ti/Si，卡在势函数上：Ti–Si、Ti–Au 只有手挑的 Lennard-Jones
交叉项，等于给界面强加了一个没有物理依据的相互作用。实测两种极端都失败：

| 交叉项取值 | 结果 |
|---|---|
| 强（ε=0.3 eV, σ=2.0 Å） | 50 ps 内界面互混放热 **2206 eV**，TBR 单调上涨 |
| 弱（Lorentz–Berthelot 混合） | 层间不粘，Si 板向真空膨胀 **15 Å** 解体 |

**结论：LJ 交叉项没有物理合理的中间值**，因为真实 Ti–Si 是硅化物化学键。

现在改用**专用机器学习势**：`wugai-whu/NEP-Cu-X-diamond-thermal-conductivity`
提供的 NEP（Cu/C/Ti/Cr/Mo/W，24k 参数），本身就是为 Cu/X/金刚石界面热输运训练的。
材料相应改为 **C(金刚石)/Ti/Cu/Ti/C**。

## 2. 已验证的基线（全部实测，非估计）

### 2.1 势函数体相标定

| 元素 | NEP | 实验 | 偏差 |
|---|---|---|---|
| C (diamond) | a₀ = 3.550 Å | 3.567 | **−0.5%** |
| Cu (fcc) | a₀ = 3.650 Å | 3.615 | +1.0% |
| Ti (hcp) | a = 2.950 Å, c/a = 1.56 | 2.951, 1.587 | −0.03% / −1.7% |

### 2.2 界面几何

Cu(111) 与 Ti(0001) 都是六方排布，最近邻 2.581 Å 与 2.950 Å，比值 = **8/7**（误差 0.03%）。
因此 **7×7 Ti 胞 ↔ 8×8 Cu 胞** 在 20.650 × 35.767 Å 内**双向精确共格**，无需人为应变。

生成器内置最近邻自检：任何原子对 <1.6 Å 直接中止构建。

### 2.3 界面热阻（三 seed）

| seed | TBR_int (m²K/W) | R² |
|---|---|---|
| 1001 | 7.292×10⁻¹⁰ | 0.981 |
| 1002 | 7.464×10⁻¹⁰ | 0.988 |
| 1003 | 7.229×10⁻¹⁰ | 0.982 |

**均值 7.33×10⁻¹⁰ m²K/W，离散 ±1.6%。**

### 2.4 稳定性（这是可信度的前提）

| 指标 | 数值 | 对照（原 LJ 体系） |
|---|---|---|
| 生产段 dPE/dt | **+0.002 eV/ps** | −1 ~ −6 eV/ps |
| 互混原子数 | **0** | 700+ |
| 平面法半-半差 | **8.9%** | 52% |

## 3. 三种估计量与交叉验算

| 方法 | 值 (m²K/W) | 说明 |
|---|---|---|
| `tbr_extrap.py`（**首选**） | **7.33×10⁻¹⁰** | 测 TBR_eff(z) for z=8..40 Å，线性拟合，截距 = 纯界面热阻 |
| `tbr_plane.py` ±30 Å | 1.71×10⁻⁹ | 含体相电阻，不适合当界面值 |
| `tbr_global_fit.py` | 2.17×10⁻⁹ | 同上，且受拟合窗口影响 |
| `postprocess.py`（窗口外推） | 7.08×10⁻¹⁰ | 与本体系外推值吻合 |

斜率的物理意义：反推的**声子**热导率约 4–5 W/mK，符合金属中"MD 只有声子、没有电子"
的预期，说明分解自洽。

## 4. 图案研究

界面图案的目标是降低 TBR。单纯在 Cu 侧刻空槽**只会减少接触面积**，因此本项目实现的是
**互锁（interlock）**：Ti 长出凸台插入 Cu 的沟槽，接触面积大于平面。

```bash
python3 tools/build_pair_hex.py --interlock 6 --groove-w 6 --groove-p 10.325 \
  --out-structure structures/IL_H6W6P10.data --out-geom structures/IL_H6W6P10.geom.json
```

参数：凸台高度 H、宽度 W、周期 P。界面处 <2.3 Å 的原子对自动删除。

## 5. 复现

```bash
# 1) 势函数（14 MB 模型已随仓库提供，也可重新下载）
bash tools/get_nep89.sh            # 可选，用于通用势对照

# 2) 结构（会打印共格误差与最近邻自检）
python3 tools/build_pair_hex.py --out-structure structures/PAIR_CuTi.data \
  --out-geom structures/PAIR_CuTi.geom.json

# 3) 门禁（期望 READY errors=0）
python3 tools/preflight.py --config config/params.json \
  --geom structures/PAIR_CuTi.geom.json --structure structures/PAIR_CuTi.data \
  --potentials potentials/potentials.in --dry-run

# 4) 建 run 并跑（约 1.4 h，24 ranks；**不要并发跑多个 NEP 作业**）
python3 tools/new_run.py --config config/params.json --scan-id PAIR_CUTI ...
MPIRUN_EXTRA="--cpu-set 0-23" bash tools/run_seed.sh workspace/run_v001 1001

# 5) 三种估计量 + 独立复核
python3 tools/tbr_extrap.py      --run-dir workspace/run_v001 --seed seed_1001
python3 tools/tbr_plane.py       --run-dir workspace/run_v001 --seed seed_1001 --z 30
python3 tools/tbr_global_fit.py  --run-dir workspace/run_v001 --seed seed_1001
python3 tools/verify_run.py      workspace/run_v001 seed_1001
```

## 6. 协议

v1.6（继承自已验证的 Si/Ti/Au 项目）：run 内退火 300→700→300 K（65 ps）→ 300 K 平衡
10 ps → 通量热身 80 ps（不计量）→ **计量 150 ps**。

## 7. 目录

```
tools/        build_pair_hex（几何+自检）、preflight（G0）、new_run、run_seed.sh、
              postprocess、tbr_plane、tbr_global_fit、tbr_extrap（首选）、
              verify_run（独立复核）、calibrate*（体相标定）
potentials/   CuTiC_NEP.txt（模型）、potentials.in、lattice.json、lock
config/       params.json（协议参数）、search_space.json
docs/         HANDOFF.md（状态与踩坑）、协议说明
reports/      evidence/（每次 run 的原始日志与结果）
```

## 8. 引用

- 势函数：NEP-Cu-X-diamond-thermal-conductivity (wugai-whu)
- 方法：NEP_CPU / GPUMD (Xu et al., MGE Advances 2025)
- 原流程：本仓库继承自 `hizuru254/Si-Ti-Au-NEMD`
