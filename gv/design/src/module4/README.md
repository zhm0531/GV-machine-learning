# GV02-04 模块4：基线优化（随机突变 + 贪心爬山）

本模块接在成员3提供的 `gv/design/src/models/` 后面，**复用模块3的参考序列、掩码、约束和目标函数**，不复制一份可能漂移的评分逻辑。

## 做什么

- `random_mutation`：从参考序列或当前最优序列附近随机产生 1 到 `--max-mutations` 个替换突变。
- `greedy_hill_climb`：从当前序列出发，枚举可变位点的单点替换，采用 best-improvement 贪心爬山。
- 目标：最大化模块3的 `target(s) = -aggregation_max(s)`，等价于最小化 AGGRESCAN a3v 的 7-aa 窗口最大聚集倾向。
- 约束：直接调用模块3 `models.constraints`。锁定位点从模块3的 `conservation_profile.json` 读取；长度、标准氨基酸、天然性质包络在搜索循环内检查。
- Pfam：搜索中不启动 HMMER；对最终各 run 的最佳候选批量使用 PF00741 + `--cut_ga` 复核。

## 参考序列/掩码口径

代码会动态读取模块3 `models.constraints` 的配置。你提供的模块3仓库当前使用：

- 参考 ID：`WP_007357596.1`
- 长度：72 aa
- 锁定位置：35
- 可变位置：37

注意：M1 文档曾写 `WP_073608715.1`，但提交的模块3可执行代码和模块2对账文件使用 `WP_007357596.1`。为保证接口一致，本模块跟随模块3当前实际配置；请团队确认后统一所有模块，不要只在模块4里单独改参考序列。

## 安装与运行

将本文件夹 `module4/` 复制到 `gv/design/src/module4/`。保持同级存在 `models/` 和 `data/`，例如：

```text
gv/design/
  data/
    conservation_profile.json
    pfam/PF00741.hmm
  src/
    models/
      objective.py
      constraints.py
      properties.py
    module4/
      __init__.py
      baselines.py
      run_baselines.py
```

从 `gv/design/src` 目录运行：

```bash
# 先做小规模冒烟测试
python module4/run_baselines.py --outdir module4_smoke --budgets 15 --seeds 42

# 正式预算（模块1的实验协议）
python module4/run_baselines.py --outdir module4_out --budgets 300,600,1200 --seeds 42,43,44,45,46
```

如果 `models/` 不在 `gv/design/src`，设置 `GVP_SRC` 指向包含 `models/` 的目录；若数据目录不同，按模块3约定设置 `GVP_DATA`。

若机器尚未安装 HMMER，算法仍会完成，`pfam_verification.csv` 会明确标记 `NOT_RUN`。正式验收建议安装 HMMER 3.x 后重跑或单独复核。可临时用 `--skip-pfam` 明确跳过最终复核。

## 输出

每个 run 的子目录（例如 `random_mutation_budget300_seed42/`）包含：

- `history.csv`：每次成功目标评估后的 best-so-far 历史，兼容成员5的中期报告脚本。
- `best.fasta`：该 run 的最佳序列。
- `reference.fasta`：参考序列。
- `summary.json`：本次运行统计。
- `config.json`：参考、预算、约束和目标定义。

根输出目录包含：

- `best.fasta`：所有 run 的最佳序列集合，供成员5约束验证/副作用分析。
- `reference.fasta`
- `summary.csv`：策略、预算、seed、评估数、目标改善、突变数和耗时。
- `pfam_verification.csv`：最终候选的 PF00741 扫描结果。
- `experiment_config.json`

可将输出目录直接传给成员5的代码3：
```bash
python "代码3" --history module4_out --verify verify_out/verification.csv --side sideeffect_out/side_effect.csv --out midterm_report.md
```

## 评估预算口径

`budget` 是**最多成功通过快速约束、实际调用模块3目标函数的次数**，包括参考序列初始评估；硬约束拒绝的候选不计入 `eval_count`。`proposal_count` 另外记录候选提议次数。贪心爬山可能因达到局部最优而在预算用完前提前停止；报告中应如实保留实际评估次数和停止原因。

每个 run 内有序列缓存，重复候选不重复调用评分器；因此评估次数代表真实新序列评分数。各策略共用同一参考序列、目标函数、约束和预算口径。

## 限制说明

- 当前目标是模块3提供的确定性理化性质代理 `properties.aggregation`，不是在线调用 AGGRESCAN 网站/API；模块3代码明确以 AGGRESCAN a3v 量表实现 `aggregation_max`。结果应称为“基于项目内 a3v 量表的计算优化”，不要声称已经完成实验验证。
- Pfam 复核仅对最终最佳候选进行，不在内循环中执行，避免秒级 HMMER 子进程拖慢搜索。
- 模块3目标函数是不可微的离散黑盒，因此这两个无梯度基线适用。


## 本次更新：对齐成员1 mod1 正式掩码

本版本以成员1 的 `mod1/outputs/readme1` 和 `mod1/outputs/mask.txt` 为准，随包附带 `config/official_mask.txt`。参考序列为 `WP_007357596.1`（72 aa）；锁定 35 个位点，可改 37 个位点。启动时会把该正式掩码同步到模块3约束层，并检查参考序列 ID/长度；如果参考序列不一致或掩码文件缺失/格式错误，会直接报错，不会静默沿用旧掩码。

正式锁定位点（0 基索引）：

```text
0, 9, 10, 11, 12, 13, 14, 16, 17, 19, 20, 21, 22, 24, 26, 32, 34, 35, 36, 37, 42, 43, 44, 46, 48, 49, 50, 52, 54, 56, 57, 58, 59, 61, 62
```

启动日志和 `experiment_config.json` 会记录锁定/可改位点数量及正式掩码来源。其余模块3目标函数、天然理化性质包络等仍复用成员3代码。
