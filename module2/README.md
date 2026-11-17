# 模块 2｜保守/可变区间 + 域/接口标注（M2 分析）交付物包

- 责任人：**成员 4**（本阶段任务调整为模块 2）
- 对应分工原文：*产出保守度图谱、可变区域清单；标注 PF00741 域边界与 GvpC 组装接口位点（供 E4 副作用审计）。交付：保守/可变区间表、域/接口位点表、可视化图。*
- 一键复现：`bash scripts/run_all.sh`（依赖已装在本目录 `.venv`）

---

## 1. 交付物对照

| 分工要求的交付物 | 本包文件 | 说明 |
|---|---|---|
| **保守/可变区间表** | `out/conservation_regions.csv` | 29 个连续区间，含列区间、参考坐标、平均/最大熵 |
| **域/接口位点表** | `out/domain_interface_sites.csv` | PF00741 域（高置信）、GvpA↔GvpC 界面候选、GvpC 重复块、往届图参考 |
| **可视化图** | `figures/fig1_conservation_curve.png`、`fig2_variability_heatmap.png`、`fig3_region_map.png` | 保守度曲线 / 一致性热图 / 区间图 |
| （分工提到的）保守度图谱 | **复用成员3 的** `gv/design/data/conservation_profile.json`（仓库 `GV02-04-member3` 分支，参考 `WP_007357596.1`）；本包提供**同参考序列**下的独立交叉验证版本 `out/conservation_profile.json` | 分工"复用 conservation_profile.json"指的就是成员3 该文件；两条独立方法已对账（见 §8） |
| **（分工要求的）可变区域清单（团队口径，37 位）** | `out/variable_regions_team_definition.csv`、`out/variable_regions_summary.json` | 37 位 = 72 − 35(成员3 掩码)；含域内外、熵标签、`majority_freq` 与口径差异标注 |
| （附加）可变位点清单（熵口径） | `out/variable_positions.csv` | 104 行（域内高变 10 + 高熵插入列 94），用于 E2 掩码敏感性 |
| 附加：掩码依据 | `out/mask_candidates_v2.json` | 基于参考坐标的 20/40/60% 秩次掩码候选（供模块 1） |
| 附加：按参考坐标的位置表 | `out/ref_position_profile.csv` | 72 行，逐位标注域内外/标签/置信度（供模块 3/5） |
| 附加：与成员3 掩码对账 | `out/reconcile_member3.json`、`figures/fig4_mask_comparison.png` | Jaccard / 精确率 / 召回率 / Spearman（见 §8） |
| 附加：域统计与 GvpC 重复统计 | `out/domain_stats.json`、`out/gvpc_repeat_summary.json` | — |

---

## 2. 方法与工具（可复现）

| 步骤 | 工具 | 命令/脚本 |
|---|---|---|
| 域 HMM | Pfam **PF00741**（Gas_vesicle, LENG=39），InterPro 下载 | `scripts/03_hmm_annotate.py` |
| 域标注 | `pyhmmer.hmmsearch --cut_ga`(=gathering)，**cpus=1 串行** | 同上 → `out/PF00741_hits.csv` |
| 主比对 | `pyhmmer.hmmalign`（HMM 引导，`trim=False`） | `scripts/04_hmmalign.py` |
| 保守度 | 逐列 Shannon 熵（忽略 gap）、gap 比例、覆盖度、滑窗 w=5/7 | `scripts/05_conservation_profile.py` |
| GvpC 重复 | Pfam **PF01304**（Gas_vesicle_C, LENG=33）；**PF05465** 作阴性对照 | `scripts/06_gvpc_interface.py` |
| 图 | matplotlib | `scripts/07_figures.py` |
| 对照比对 | FAMSA（`pyfamsa`），**仅作对照**，见下"局限" | `scripts/02_build_msa.py` |

数据来源：`pan_resources/gv/design/data/GvpA.fasta`（856 条）、`GvpC.fasta`（3225 条）。

---

## 3. 关键结果

### 3.1 GvpA（参考序列 `WP_007357596.1`，72 aa —— 与成员3 的 profile 同一参考）

| 指标 | 结果 |
|---|---|
| PF00741 域命中率 | **854/856**（gathering 阈值） |
| 参考序列域区间 | **11–49**（39 aa，E=1.03e-18，score=62.9） |
| 主比对 | 856 × 153 列（**0 个全 gap 列**；39 个 HMM 匹配列） |
| 匹配列熵 | min 0.00 / 中位 **0.31** / max 1.68；**严格保守列 2 个** |
| 域内保守 / 中间 / 可变 | **10 / 19 / 10** |
| 侧翼（域外） | 33 位，高熵（N 端 1–10、C 端 50–72 均为超变区） |
| 按参考坐标的可变位点 | 域内高变 10 + 侧翼 33 = **43 位** |

**对分工文档"可变区域约 37 位"的核对（已找到准确来源）**：37 = **72 − 35**，即成员3 的 `conservation_profile.json` 用 `majority_freq ≥ 0.9` 锁定了 **35 个保守位**，剩余 **37 位**即分工所说的"可变区域"。这**不是**往届伪比对的产物（先前版本的本 README 曾误判为后者，已更正）。本包按熵口径得到 43 位（10 个域内高变 + 33 个侧翼），差异来自"侧翼是否计入可变区"的定义，报告里应写明采用哪个口径。

### 3.2 GvpC（3225 条）

| 指标 | 结果 |
|---|---|
| PF01304 重复域命中 | **179/3225（5.5%）**，共 884 个域 |
| 每序列重复单元 | 中位 **5 个 × 32 aa** |
| 重复块 | 占序列中位 **84%**，起点中位 19，块长中位 150 aa |
| PF05465（古菌） | 0 命中（阴性对照符合预期） |

### 3.3 两条重要更正（影响全组）

1. **口径更正**：项目池写的是 `pfam01132`，实为 **Elongation factor P (EF-P) OB domain**，与气囊蛋白无关；正确条目是 **PF00741**（Gas vesicle protein）。本包已按 PF00741 执行，建议全组统一并在报告注明。
2. **往届 `.a3m` 不是比对**：`create_simple_msa.py` 只是"截断到最短序列长度"（856 条 × 60 列、零 gap、全部等于原序列前 60 残基）。本包改用 HMM 引导比对重建，结论更可靠。

---

## 4. 字段说明（下游消费用）

```
conservation_profile.json
  reference{id,length} / alignment{n_seq,n_cols,n_match_cols} / method
  thresholds{match_col_entropy_quantiles{q10,q25,q50,q75,q90}, conserved, variable}
  windows{w5:[...],w7:[...]}
  columns[{col, ref_pos, col_type(match|insert), consensus, entropy,
           gap_frac, coverage, confidence(high|medium|low), label}]

conservation_regions.csv
  region_id, label(conserved|intermediate|variable|insert),
  start_col, end_col, length, ref_start, ref_end, mean_entropy, max_entropy

domain_interface_sites.csv
  site_id, protein, start, end,
  type(PF00741_domain|assembly_interface|prior_work_reference),
  evidence, confidence(high|medium|low), source

ref_position_profile.csv
  ref_pos, col, col_type, consensus, entropy, label, in_PF00741_domain, confidence
```

---

## 5. 与其它模块的接口

| 模块 | 从本包取什么 |
|---|---|
| 模块 1（正式掩码） | `mask_candidates_v2.json`（20/40/60% 秩次掩码）、`ref_position_profile.csv`（逐位熵与标签）——**注意：GvpA 无"熵=0 即锁定"的位点，必须用秩次法或 PF00741 域内位点** |
| 模块 3（约束层） | `domain_interface_sites.csv`（域区间作硬约束/包络）、`ref_position_profile.csv`（两级判定的位点清单） |
| 模块 5（E4 审计） | `domain_interface_sites.csv`（接口位点改动数）、`conservation_regions.csv`（意外改动是否落在可变区） |

---

## 6. 局限与风险（必须在报告中如实说明）

1. **GvpC 侧界面为推定**：PF01304 只覆盖 179/3225 条序列；且**无实验复合物结构**，GvpA 侧的残基级接触面未确定（`confidence=low`）。往届 region 图（`reference/`）未复核，仅作参考。
2. **插入列（侧翼）的比对是"打包式"**：HMM 引导比对把长度可变的 N/C 端放入插入状态，其列熵可反映变异性，但不能当作严格的位置同源。
3. **FAMSA 作为对照失败**：191 列中 93 列为 100% gap（对短而分歧的家族过比对），故未采用；已在 `out/GvpA_famsa.aln` 保留供核查。
4. **域边界为家族级**：PF00741 是 family 条目（M=39），对 60–100 aa 的 GvpA 只覆盖约 46%；域外序列并非"无功能"，只是无该家族模型覆盖。
5. **HTTPS/网络依赖**：HMM 文件来自 InterPro，若离线需预先放入 `pfam/`。

---

## 7. 目录

```
module2/
├── README.md                  ← 本文件
├── pfam/                      PF00741.hmm / PF01304.hmm / PF05465.hmm
├── scripts/                   02~09 管线 + run_all.sh（+ fetch_pan.py 取数）
├── out/                       全部数据交付物（json/csv）
├── figures/                   4 张可视化图
└── reference/                 网盘参考图（GvpA/GvpC region 五联图等）+ gvpc_metadata.csv
```

---

## 8. 与 `GV02-04-member3` 分支的关系与对账（重要）

**结论：不是冲突，而是互相印证。** 分工里模块2 要求"复用 conservation_profile.json"，指的就是成员3 在 `GV02-04-member3` 分支产出的
`gv/design/data/conservation_profile.json`（参考 `WP_007357596.1`，`mask_threshold=0.9`，锁定 35 位）。
本包已**改用同一参考序列**重跑全部交付物，并与其掩码逐位对账（`scripts/08_reconcile_member3.py`）：

| 对账项 | 结果 | 解读 |
|---|---|---|
| 成员3 掩码 | 35 锁定 / 37 自由（majority_freq ≥ 0.9） | 37 = 分工"约 37 位"的来源 |
| 本包 top40（29 位） vs 成员3 | 交集 **29**，精确率 **1.000**，召回率 0.829，Jaccard 0.829 | **本包 top40 是成员3 掩码的完全子集** → 两条独立方法一致 |
| 本包 top20（14 位） vs 成员3 | 精确率 **1.000**，召回率 0.40 | 更严格的核心里全部落在成员3 掩码内 |
| 本包 top60（43 位） vs 成员3 | 精确率 0.744，召回率 0.914 | 放宽到 60% 时会多锁成员3 未锁的位点 |
| majority_freq vs 熵（取负） | **Spearman = 0.930** | 两种保守性度量高度一致 |
| 成员3 锁定位分布 | 域内（11–49）**23** 位，域外 **12** 位 | 分歧集中在侧翼：成员3 按频率锁、本包按秩次排 |

**给团队的建议**：
1. **掩码以成员3 的 35 位为准**（模块3 的 `constraints.py` 已在消费它），本包 top40 作为其子集可用于"严格掩码"敏感性对照；
2. 侧翼那 12 位争议位点建议在报告中列为"口径差异"，并用 E2 的掩码敏感性实验（top20/40/60）给出实证；
3. 另一处仍需注意：成员3 分支目录里同时放了 `PF01132.hmm`，而 **PF01132 = Elongation factor P (EF-P) OB domain，与 GvpA 无关**，不应使用；
   GvpA 的正确条目是 **PF00741**（本包已按此执行）。

**两条独立路径的一致结论**：包内 `GvpA_msa.a3m` 是 `create_simple_msa.py` 截断到最短长度(60)的伪比对，不能用于保守性统计
（成员3 的脚本注释中也独立记录了同一结论）。

**重锚定方式**：脚本已支持环境变量，无需改代码——
```bash
GVP_REF_ID=WP_007357596.1 bash scripts/run_all.sh   # 团队参考（默认）
GVP_REF_ID=WP_073608715.1 bash scripts/run_all.sh   # 换回另一条参考做敏感性
```
