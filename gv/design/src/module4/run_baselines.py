# -*- coding: utf-8 -*-
"""命令行入口：运行随机突变与贪心爬山的预算/随机种子对比实验。"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_SRC = HERE.parent
for p in (os.environ.get("GVP_SRC", ""), str(DEFAULT_SRC)):
    if p and os.path.isdir(os.path.join(p, "models")) and p not in sys.path:
        sys.path.insert(0, p)

from baselines import run_baseline, REF_ID, REF_SEQ, MASK, FREE


def _write_fasta(path, records):
    with open(path, "w", encoding="utf-8") as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")


def _write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def _verify_pfam(results, outdir, skip=False):
    """对最终每次 run 的 best 做 PF00741 批量复核；无 HMMER 时明确标记未执行。"""
    path = os.path.join(outdir, "pfam_verification.csv")
    rows = []
    if skip:
        rows = [{"name": f"{r.strategy}_b{r.budget}_s{r.seed}",
                 "passed": "SKIPPED_BY_USER", "evalue": "", "score": "",
                 "status": "用户指定跳过"} for r in results]
    else:
        try:
            from models.constraints import pfam_scan
            names = [f"{r.strategy}_b{r.budget}_s{r.seed}" for r in results]
            hits = pfam_scan([r.best_sequence for r in results], names=names)
            for h in hits:
                rows.append({
                    "name": h["name"], "passed": h["passed"],
                    "evalue": h["evalue"], "score": h["score"],
                    "status": "executed",
                })
        except Exception as e:
            rows = [{"name": f"{r.strategy}_b{r.budget}_s{r.seed}",
                     "passed": "NOT_RUN", "evalue": "", "score": "",
                     "status": f"未执行：{type(e).__name__}: {e}"} for r in results]
            print(f"[警告] PF00741 复核未执行：{type(e).__name__}: {e}", file=sys.stderr)
    _write_csv(path, rows)
    return path


def main():
    ap = argparse.ArgumentParser(description="GV02-04 模块4：随机突变 + 贪心爬山")
    ap.add_argument("--outdir", default="module4_out", help="结果输出目录")
    ap.add_argument("--budgets", default="300,600,1200", help="评估预算，逗号分隔")
    ap.add_argument("--seeds", default="42,43,44,45,46", help="随机种子，逗号分隔")
    ap.add_argument("--strategies", default="random_mutation,greedy_hill_climb",
                    help="可选 random_mutation,greedy_hill_climb")
    ap.add_argument("--max-mutations", type=int, default=3,
                    help="随机突变每次改变的可变位点数上限")
    ap.add_argument("--skip-pfam", action="store_true",
                    help="跳过最终 PF00741/HMMER 复核（默认尝试执行）")
    args = ap.parse_args()

    budgets = [int(x.strip()) for x in args.budgets.split(",") if x.strip()]
    seeds = [int(x.strip()) for x in args.seeds.split(",") if x.strip()]
    strategies = [x.strip() for x in args.strategies.split(",") if x.strip()]
    allowed = {"random_mutation", "greedy_hill_climb"}
    if not budgets or any(b < 1 for b in budgets):
        ap.error("--budgets 必须是正整数列表")
    if not seeds:
        ap.error("--seeds 不能为空")
    if not strategies or any(s not in allowed for s in strategies):
        ap.error(f"--strategies 只能包含 {sorted(allowed)}")
    if args.max_mutations < 1:
        ap.error("--max-mutations 必须 >= 1")

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    print(f"参考序列：{REF_ID}，长度 {len(REF_SEQ)} aa")
    print(f"锁定位点：{len(MASK)}；允许突变位点：{len(FREE)}")
    print("目标：maximize -aggregation_max（等价于最小化 7-aa 窗口聚集热点）")
    print(f"实验：strategies={strategies}, budgets={budgets}, seeds={seeds}")

    results = []
    for strategy in strategies:
        for budget in budgets:
            for seed in seeds:
                run_name = f"{strategy}_budget{budget}_seed{seed}"
                run_dir = os.path.join(outdir, run_name)
                print(f"\n[{len(results)+1}] {run_name} ...", flush=True)
                kwargs = {"max_mutations": args.max_mutations} if strategy == "random_mutation" else {}
                result = run_baseline(strategy, budget, seed, outdir=run_dir, **kwargs)
                results.append(result)
                delta = result.best_score - result.reference_score
                print(
                    f"  eval={result.eval_count}/{budget}, proposals={result.proposal_count}, "
                    f"invalid={result.invalid_count}, edits={sum(a!=b for a,b in zip(result.best_sequence, REF_SEQ))}, "
                    f"aggregation_max={-result.reference_score:.4f}->{-result.best_score:.4f}, "
                    f"Δtarget={delta:+.4f}, stop={result.stop_reason}, "
                    f"time={result.elapsed_seconds:.2f}s"
                )

    # 模块5的代码3可直接读取该根目录（递归搜索各 run 子目录中的 history.csv）。
    # 另写 summary.csv 方便表格汇总；不在根目录生成同名 history.csv，避免重复读取。
    summary_rows = []
    all_best = []
    for r in results:
        summary_rows.append({
            "strategy": r.strategy, "budget": r.budget, "seed": r.seed,
            "ref_id": r.ref_id, "eval_count": r.eval_count,
            "proposal_count": r.proposal_count, "invalid_count": r.invalid_count,
            "cache_hits": r.cache_hits, "reference_score": r.reference_score,
            "best_score": r.best_score, "target_improvement": r.best_score-r.reference_score,
            "reference_aggregation_max": -r.reference_score,
            "best_aggregation_max": -r.best_score,
            "n_edits": sum(a != b for a, b in zip(r.best_sequence, REF_SEQ)),
            "elapsed_seconds": round(r.elapsed_seconds, 4), "stop_reason": r.stop_reason,
        })
        all_best.append((f"{r.strategy}|budget={r.budget}|seed={r.seed}", r.best_sequence))
    _write_csv(os.path.join(outdir, "summary.csv"), summary_rows)
    _write_fasta(os.path.join(outdir, "best.fasta"), all_best)
    _write_fasta(os.path.join(outdir, "reference.fasta"), [(REF_ID, REF_SEQ)])

    verify_path = _verify_pfam(results, outdir, skip=args.skip_pfam)
    with open(os.path.join(outdir, "experiment_config.json"), "w", encoding="utf-8") as f:
        json.dump({
            "project": "GV02-04",
            "module": 4,
            "reference_id": REF_ID,
            "reference_length": len(REF_SEQ),
            "locked_positions_count": len(MASK),
            "free_positions_count": len(FREE),
            "strategies": strategies, "budgets": budgets, "seeds": seeds,
            "max_mutations": args.max_mutations,
            "objective": "maximize -aggregation_max",
            "constraint_source": "module1 mod1/outputs/mask.txt (bundled as module4/config/official_mask.txt), with module3 constraints/objective",
            "official_mask_zero_based": list(MASK),
            "pfam_verification_csv": os.path.basename(verify_path),
            "notes": [
                "The formal locked-site mask from member1 mod1 is authoritative; module3 MASK/FREE are synchronized at import.",
                "Official mask: 35 locked sites and 37 mutable sites; only FREE sites are mutated.",
                "Infeasible candidates do not consume objective evaluation budget.",
                "PF00741 HMMER scan is performed only on final best candidates.",
                "Reference sequence follows mod1 readme1: WP_007357596.1, length 72 aa.",
            ],
        }, f, ensure_ascii=False, indent=2)

    print(f"\n完成。输出目录：{outdir}")
    print(f"汇总：{os.path.join(outdir, 'summary.csv')}")
    print(f"最佳序列：{os.path.join(outdir, 'best.fasta')}")
    print(f"PF00741 复核：{verify_path}")


if __name__ == "__main__":
    main()
