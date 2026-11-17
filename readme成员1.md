# 成员5 · 模块5 说明（GV02-04 定向序列优化）

本模块负责 **M3 收尾：验证 + 副作用分析 + 中期演示**。各成员（模块1–4）产出的优化序列交由本模块验收：是否合规、有没有引入副作用、能否汇总成中期材料。

## 所属分工与实现

| 职责 | 对应阶段 | 实现 |
|---|---|---|
| 约束验证（合规性） | 模块5 – 块1 | `module5/verify.py` |
| 副作用分析（E4） | 模块5 – 块2 | `module5/side_effect.py` |
| 中期演示材料生成 | 模块5 – 块3 | `module5/report.py` |
| 一键串起块1–3 | 模块5 | `run_all.py` |

> 说明：核心逻辑统一收敛到 `module5/` 包内；共享工具（fasta/掩码/数据路径/性质加载）集中在 `module5/common.py`，避免重复实现。

## 目录结构

```
jiqixuexizongheshijian/
├─ run_all.py           # 一键运行全流程
├─ readme成员5.md
└─ module5/                # 核心实现（Python 包）
   ├─ __init__.py
   ├─ common.py         # 共享：fasta/掩码/数据路径/性质加载
   ├─ verify.py         # 块1 约束验证
   ├─ side_effect.py    # 块2 副作用分析
   └─ report.py         # 块3 中期报告
```

## 输入与输出

### 输入（依赖模块1–4 产出）
| 文件 | 说明 | 来源 |
|---|---|---|
| `best.fasta` | 优化器产出的候选序列 | 模块4 |
| `mask.txt` / `conservation_profile.json` | 锁定位点 / 保守性图谱 | 模块1 |
| `history.csv` | 优化收敛历史（可给目录） | 模块4 |
| `GVP_DATA` / `--data` 指向的 gv 数据目录 | 数据与性质模块 | 项目 |

### 输出
| 模块 | 输出 | 含义 |
|---|---|---|
| `module5/verify.py` | `verify_out/verification.csv`、`verify_out/constraint_report.txt` | 约束验证明细/报告 |
| `module5/side_effect.py` | `sideeffect_out/side_effect.csv`、`sideeffect_out/side_effect_report.txt` | 副作用分析明细/报告 |
| `module5/report.py` | `midterm_report.md`（可选 `convergence.png`） | 中期报告/演示素材 |

## 运行方法

在 `jiqixuexizongheshijian` 目录下运行：

```powershell
# 一键全流程
python run_all.py --input best.fasta --data D:/aineeding/gv/design/data

# 分步运行
python -m module5.verify --input best.fasta --ref ref.txt --mask mask.txt --outdir verify_out
python -m module5.side_effect --input best.fasta --data D:/aineeding/gv/design/data --mask mask.txt
python -m module5.report --verify verify_out/verification.csv --side sideeffect_out/side_effect.csv --out midterm_report.md
```

> 数据目录也可用环境变量 `GVP_DATA` 指定；源码目录用 `GVP_SRC` 指定。

## 关键参数（改动前先确认）
| 参数 | 值 | 说明 |
|---|---|---|
| `VARIABLE_THR` | 0.70 | 多数频率≤它视为可变位点 |
| `OUTSIDE_RATIO_LIMIT` | 0.10 | 非可变区改动占比上限 |
| `FOLD_DROP_LIMIT` | 0.05 | 折叠代理下降上限 |
| 参考序列 | `WP_007357596.1`（72aa，锁定35/可改37） | 以 `conservation_profile.json` 为准 |

## 待办事项
- [ ] 模块4 产出 `best.fasta` 与 `history.csv` 后，把真实数据传进来重跑。
- [ ] 参考序列口径以 `conservation_profile.json` 的 `ref_id`（`WP_007357596.1`）为准。
- [ ] 若装了 HMMER / ColabFold，替换重判定 / 折叠代理以启用完整验证。
- [ ] 汇总全体模块，产出最终中期检查演示材料。

## 说明
- 全部**纯标准库**实现，无第三方依赖（`matplotlib` 仅 `module5/report.py` 的可选出图）。
- HMMER、AlphaFold 未安装时相关功能**优雅降级**而非报错。
- `module5/verify.py` / `side_effect.py` / `report.py` 均可独立运行（`python -m module5.<模块>`），也可由 `run_all.py` 一键串起。
- 早期曾以 `代码1/2/3` 三个独立脚本提供入口，现已并入 `module5/` 包，命令行为 `python -m module5.*`。