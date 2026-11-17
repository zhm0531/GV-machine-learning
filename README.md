# GV02-04 定向序列优化

《机器学习综合实践》Track 1（GV，气囊蛋白）T02 项目。

给定蓝藻气囊主壳体蛋白 **GvpA** 参考序列，在**保持 Pfam 域约束与关键位点保守**的
前提下，定向优化目标性质，并比较不同优化策略的效果。

三个研究问题（课程统一框架，见 [`GV02-04-项目要求整合.md`](GV02-04-项目要求整合.md)）：

* **RQ1** 保持 Pfam 域约束满足的前提下，性质提升的极限在哪里？提升幅度与序列变化量的关系？
* **RQ2** 不同优化策略（梯度引导 / 强化学习 / 进化策略 / 贝叶斯优化）在收敛速度、样本效率、约束满足率上如何对比？
* **RQ3** 优化后的序列是否引入新问题（功能变化、结构不稳定、组装接口不兼容）？

---

## 目录结构

```
gv/design/
  data/
    GvpA.fasta                  856 条天然 GvpA 序列
    conservation_profile.json   保守位点图谱 + 天然性质包络（生成物，已入库）
    pfam/PF00741.hmm            Pfam 官方模型
  preprocessing/
    build_conservation_profile.py    生成图谱与包络
  src/
    models/
      properties.py           理化性质泛函（纯标准库）
      constraints.py          硬约束层：Pfam + 锁定 + 天然包络
      objective.py            目标函数 + 爬山基线
    score/                    官方评分模块（依赖外部工具，见「注意事项」）
    experiments/              对比实验（待填充）
```

---

## 环境依赖

* **Python 3.10+**，仅用标准库，无需第三方 Python 包
* **HMMER 3.x**（提供 `hmmsearch`），用于 Pfam 约束校验

```bash
# Debian / Ubuntu
sudo apt-get install hmmer
# macOS
brew install hmmer
```

可选（仅 M4 结构类副作用分析需要，非本项目运行前提）：
ColabFold / AlphaFold、PyRosetta、MAFFT、MMseqs2。

---

## 快速开始

```bash
cd gv/design

# 1) 生成保守位点图谱与天然性质包络
python3 preprocessing/build_conservation_profile.py

# 2) 分层自检
cd src
python3 models/properties.py     # 理化性质
python3 models/constraints.py    # 硬约束（含 Pfam 扫描）
python3 models/objective.py      # 目标函数 + 爬山基线端到端
```

作为库调用：

```python
import sys; sys.path.insert(0, "gv/design/src")
from models.objective import score, feasible_score, objective, verify
from models.constraints import REF_SEQ, MASK, FREE, ConstraintViolation

score(REF_SEQ)          # 目标 + 全约束；违规抛 ConstraintViolation
feasible_score(seq)     # 不可行返回 -inf，供优化器直接调用
objective(seq)          # 目标 + 全部辅助指标
verify([seq], ["x"])    # Pfam00741 官方模型批量复核
```

关键路径可用环境变量覆盖：

| 变量 | 默认 | 含义 |
|---|---|---|
| `GVP_DATA` | `gv/design/data` | 数据目录 |
| `GVP_REF_ID` | `WP_007357596.1` | 参考序列 ID |
| `GVP_MASK_THR` | `0.9` | 锁定位点的多数残基频率阈值 |
| `GVP_ENVELOPE_MODE` | `minmax` | 包络口径，可选 `p5p95` |
| `GVP_ENVELOPE_SKIP` | `aggregation_max` | 从包络约束中剔除的性质（即目标项） |

---

## 目标函数与约束

$$
\max_{s}\ f(s) \quad \text{s.t.} \quad
s \in \mathcal{C}_{\text{struct}} \cap \mathcal{C}_{\text{pfam}} \cap \mathcal{C}_{\text{env}}
$$

**约束是前提条件，不是目标函数中的加权项。** 违反约束的候选直接判为不可行、
不参与排序，而不是扣分。

| 层 | 内容 | 实现 |
|---|---|---|
| 目标 $f$ | $-\phi_{\text{max}}$，即**聚集热点最小化**（AGGRESCAN a3v，7 残基窗最大值） | `properties.aggregation` |
| 约束 1 | Pfam00741 + gathering threshold（官方 profile HMM） | `constraints.pfam_scan` |
| 约束 2 | 35 个锁定位点不可改（位点多数残基频率 ≥ 0.9） | `constraints.check` |
| 约束 3 | 5 项理化性质落在 856 条天然序列实测范围内 | `constraints.check_properties` |
| 辅助指标 | conserv / novelty / μH / β倾向 / GRAVY / pI / 净电荷 | 报告与 RQ3 用 |

天然性质包络（由 856 条天然序列实测，无主观阈值）：

| 性质 | 下界 | 上界 |
|---|---|---|
| `aggregation_mean` | 0.2441 | 0.3627 |
| `gravy` | 0.4125 | 0.8667 |
| `beta_propensity` | 1.0271 | 1.0987 |
| `pI` | 3.9941 | 6.5075 |
| `net_charge_pH7` | −7.7717 | −0.7770 |

`aggregation_max` 不在包络内 —— 它正是被优化的目标项。若既约束又优化同一个量，
最优解恒等于该约束边界，RQ1 的「提升极限」会退化成恒等式。

选聚集倾向作为目标的依据：野生型 $\phi_{\text{max}} = 0.7134$，天然序列最低
0.6764，优化方向明确且提升空间真实；而 GRAVY 距上界仅 0.09，两亲性类目标基本无空间。

理论推导（含三个被否定设计的失败机理）见
[`docs/目标函数设计与理论模型.md`](docs/目标函数设计与理论模型.md)。

---

## 实测结果

```
野生型   聚集热点 0.7134   conserv 0.6141   GRAVY 0.7750   pI 4.37
爬山 11 步（净改动 6 位：23 I->D, 28 W->C, 33 L->A, 38 L->D, 45 I->C, 60 V->A）
聚集热点 0.7134 -> 0.5073   (-28.9%)   天然最低 0.6764
conserv  0.6141 -> 0.3837
Pfam00741  野生型 E=2.40e-21 passed=True
           优化后 E=2.00e-14 passed=True
距包络边界余量  aggregation_mean=0.0092  beta_propensity=0.0080
                gravy=0.1083  pI=0.1359  net_charge_pH7=2.6957
```

* **RQ1 初步回答**：提升极限由**交织的包络约束**决定 —— 收敛点处
  `beta_propensity` 与 `aggregation_mean` 余量仅约 0.008–0.009，已接近绑定。
* **RQ3 现成线索**：优化后 `conserv` 由 0.6141 降至 0.3837，序列明显偏离天然家族，
  需在副作用分析中量化。

---

## 注意事项

### 1. 任务书中的 Pfam 编号有误

任务书撰写了「保守域 pfam01132」。实测：

| 编号 | 真实名称 | 扫 856 条 GvpA |
|---|---|---|
| PF01132 | Elongation factor P (EF-P) OB domain | **0 命中** |
| **PF00741** | **Gas vesicle protein**（含 GvpA / GvpJ） | 854/856 通过，参考序列 E=2.4e-21 |

**本项目一律使用 PF00741。**

> E-value 随搜索空间大小变化：同一序列单条搜索得 1.2e-21，在 856 条库中搜索得
> 2.4e-21。此处统一报告库搜索值。

### 2. `data/GvpA_msa.a3m` 不是比对结果

该文件的 856 条记录长度全部为 60、gap 数 0、小写残基数 0 —— 由
`src/create_simple_msa.py` 将每条序列**截断到最短长度**拼接而成，不含任何比对信息。
用它按列统计保守性得到的是截断伪影。

因此本项目用 `preprocessing/build_conservation_profile.py` 重建：以参考序列为锚，
对 856 条序列做 Needleman-Wunsch 全局比对后堆叠成位点图谱。纯标准库实现，
结果可重复（重建逐字节一致）。

> 若需含系统发生信息的真 MSA，需安装 MAFFT / ClustalOmega 另行构建。

### 3. 官方 `src/score/` 模块不可直接作为优化目标

`score/score_net.py` 在缺少 ColabFold / MMseqs2 / PyRosetta 时会整体回退到
`score_mock`，其中含 `random.uniform` —— 同一条序列两次调用会得到不同分数
（实测 0.6835 与 0.8050），无法用于优化。`score/score_conserv.py` 另需
`Bio.AlignIO`，且读取的是上述非比对文件。

本项目的 `models/` 与之独立，不依赖外部工具。

---

## 数据来源

课程网盘 `resources` 目录（GV / PR / XAI 共用），含 856 条天然 GvpA 序列等原始数据。

> 该网盘为课程内部资源，公开仓库中已隐去访问链接与提取码；需要原始数据的读者
> 可通过 Pfam / UniProt 自行检索 GvpA 家族，或联系仓库作者说明用途。

Pfam 模型 `PF00741.hmm` 取自 InterPro API：

```
https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00741/?annotation=hmm
```

---

## 相关文档

| 文档 | 内容 |
|---|---|
| [`GV02-04-项目要求整合.md`](GV02-04-项目要求整合.md) | 课程任务书原文整合（含 RQ、里程碑、评分标准） |
| [`docs/目标函数设计与理论模型.md`](docs/目标函数设计与理论模型.md) | 目标函数的形式化推导、三个失败设计的机理、可解性分析 |

---

## 许可

本项目代码以 [MIT License](LICENSE) 发布。数据集与第三方模型（Pfam `PF00741.hmm` 等）
版权归各自来源所有，遵循其原始许可。

