# -*- coding: utf-8 -*-
"""GV02-04 模块4：基线优化算法。

直接复用模块3的 models.objective / models.constraints：
- 最大化 target(s) = -aggregation_max(s)，因此分数越高越好；
- 只在 FREE（非锁定）位点生成突变；
- 结构/天然性质硬约束由模块3的 constraints.check 负责；
- Pfam PF00741 的 HMMER 扫描只在最终候选上批量复核，不进入搜索内循环。

本文件不使用随机 mock 评分。每个合法、未缓存的新候选只调用一次真实的
模块3目标函数，并计入 eval_count。不可行候选不会计入目标评估预算。
"""
from __future__ import annotations

import csv
import json
import os
import random
import time
from dataclasses import dataclass, asdict
from typing import Callable, Dict, List, Optional, Tuple
from pathlib import Path

# run_baselines.py 会将 gv/design/src 放入 sys.path；也支持从该目录直接 import。
try:
    # 先加载模块3约束模块，再以成员1 mod1 的正式掩码作为唯一权威配置。
    from models import constraints as constraints_module
except ImportError as exc:
    raise ImportError(
        "无法导入模块3。请从 gv/design/src 运行，或设置 GVP_SRC 指向包含 models/ 的目录。"
    ) from exc

_MASK_FILE = Path(__file__).resolve().parent / "config" / "official_mask.txt"
def _read_official_mask(path: Path) -> Tuple[int, ...]:
    if not path.exists():
        raise FileNotFoundError(f"缺少模块1正式掩码文件：{path}")
    values = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        values.append(int(line))
    mask = tuple(values)
    if len(mask) != 35 or len(set(mask)) != 35 or tuple(sorted(mask)) != mask:
        raise ValueError("正式掩码格式错误：应为 35 个不重复且升序的 0 基索引")
    if mask[0] < 0 or mask[-1] >= len(constraints_module.REF_SEQ):
        raise ValueError("正式掩码索引超出参考序列范围")
    return mask

OFFICIAL_MASK = _read_official_mask(_MASK_FILE)
if constraints_module.REF_ID != "WP_007357596.1" or len(constraints_module.REF_SEQ) != 72:
    raise RuntimeError(
        f"参考序列与模块1正式掩码不一致：当前 {constraints_module.REF_ID}/"
        f"{len(constraints_module.REF_SEQ)} aa，期望 WP_007357596.1/72 aa"
    )
# 覆盖模块3配置中的 MASK/FREE，并同步约束门面；check/satisfies 会读取模块全局变量。
constraints_module.MASK = OFFICIAL_MASK
constraints_module.FREE = tuple(i for i in range(len(constraints_module.REF_SEQ))
                                if i not in set(OFFICIAL_MASK))
constraints_module.default_constraints.mask = constraints_module.MASK
constraints_module.default_constraints.free = constraints_module.FREE

# 必须在同步正式掩码后导入 objective。
from models import objective as objective_module
objective_module.MASK = constraints_module.MASK
objective_module.FREE = constraints_module.FREE
ALPHABET = constraints_module.ALPHABET
REF_ID = constraints_module.REF_ID
REF_SEQ = constraints_module.REF_SEQ
MASK = constraints_module.MASK
FREE = constraints_module.FREE
check = constraints_module.check
satisfies = constraints_module.satisfies
default_constraints = constraints_module.default_constraints


@dataclass
class Evaluation:
    sequence: str
    score: float                 # target = -aggregation_max，越大越好
    aggregation_max: float
    aggregation_mean: float
    n_edits: int
    conserv: float
    novelty: float


@dataclass
class RunResult:
    strategy: str
    budget: int
    seed: int
    ref_id: str
    reference_sequence: str
    best_sequence: str
    best_score: float
    reference_score: float
    eval_count: int
    proposal_count: int
    invalid_count: int
    cache_hits: int
    elapsed_seconds: float
    stop_reason: str
    history: List[dict]


class Evaluator:
    """对模块3目标函数的轻量包装，按 run 统计真实可行候选评估次数。"""

    def __init__(self):
        self.obj = objective_module.GvpObjective()
        self.cache: Dict[str, Evaluation] = {}
        self.proposal_count = 0
        self.invalid_count = 0
        self.cache_hits = 0

    @property
    def eval_count(self) -> int:
        return self.obj.eval_count

    def evaluate(self, seq: str) -> Optional[Evaluation]:
        self.proposal_count += 1
        if seq in self.cache:
            self.cache_hits += 1
            return self.cache[seq]
        # 快速约束先行；违规候选不调用目标函数，也不计入真实 F 评估预算。
        if not satisfies(seq):
            self.invalid_count += 1
            return None
        details = self.obj.evaluate(seq)  # 模块3真实目标/性质计算；计数一次
        ev = Evaluation(
            sequence=seq,
            score=float(details["target"]),
            aggregation_max=float(details["aggregation_max"]),
            aggregation_mean=float(details["aggregation_mean"]),
            n_edits=int(details["n_edits"]),
            conserv=float(details["conserv"]),
            novelty=float(details["novelty"]),
        )
        self.cache[seq] = ev
        return ev


def _mutate(sequence: str, rng: random.Random, max_mutations: int = 3) -> str:
    """只在模块3的 FREE 位点上进行替换突变；不插入、不删除。"""
    if not FREE:
        return sequence
    n = rng.randint(1, min(max_mutations, len(FREE)))
    positions = rng.sample(list(FREE), n)
    chars = list(sequence)
    for pos in positions:
        choices = [aa for aa in ALPHABET if aa != chars[pos]]
        chars[pos] = rng.choice(choices)
    return "".join(chars)


def _history_row(strategy, budget, seed, ev_count, proposal_count, current, best, note=""):
    return {
        "strategy": strategy,
        "budget": budget,
        "seed": seed,
        "eval": ev_count,
        "proposal": proposal_count,
        "current_score": current.score if current else "",
        "best_so_far": best.score,
        "best_aggregation_max": best.aggregation_max,
        "best_n_edits": best.n_edits,
        "note": note,
    }


def random_mutation_search(
    budget: int = 300,
    seed: int = 42,
    max_mutations: int = 3,
    max_proposal_factor: int = 100,
) -> RunResult:
    """随机突变基线：随机生成 1..max_mutations 个可变位点的候选。

    budget 是最多调用多少次真实目标评分（包括参考序列的初始评估）。
    无效候选不计入该预算；为防止一直采到不可行候选，另设 proposal 上限。
    """
    if budget < 1:
        raise ValueError("budget 必须 >= 1")
    if max_mutations < 1:
        raise ValueError("max_mutations 必须 >= 1")
    rng = random.Random(seed)
    evaluator = Evaluator()
    t0 = time.perf_counter()

    ref = evaluator.evaluate(REF_SEQ)
    assert ref is not None, "模块3参考序列未通过自身约束，请先修复模块3数据/配置"
    current = best = ref
    history = [_history_row("random_mutation", budget, seed, evaluator.eval_count,
                            evaluator.proposal_count, current, best, "reference")]
    seen = {REF_SEQ}
    max_proposals = max(budget * max_proposal_factor, budget + 10)
    stop_reason = "evaluation_budget"

    while evaluator.eval_count < budget and evaluator.proposal_count < max_proposals:
        cand = _mutate(current.sequence if rng.random() < 0.15 else REF_SEQ,
                       rng, max_mutations)
        if cand in seen:
            # 再尝试从当前最优附近生成；若仍重复，也交给 Evaluator 记一次 proposal
            #（命中缓存时不会重复调用评分器），以保证 proposal 上限一定生效。
            cand = _mutate(best.sequence, rng, max_mutations)
        seen.add(cand)
        ev = evaluator.evaluate(cand)
        if ev is None:
            continue
        current = ev
        if ev.score > best.score:
            best = ev
        history.append(_history_row(
            "random_mutation", budget, seed, evaluator.eval_count,
            evaluator.proposal_count, current, best))
    else:
        if evaluator.proposal_count >= max_proposals and evaluator.eval_count < budget:
            stop_reason = "proposal_limit"
        else:
            stop_reason = "evaluation_budget"

    return RunResult(
        strategy="random_mutation", budget=budget, seed=seed, ref_id=REF_ID,
        reference_sequence=REF_SEQ, best_sequence=best.sequence,
        best_score=best.score, reference_score=ref.score,
        eval_count=evaluator.eval_count, proposal_count=evaluator.proposal_count,
        invalid_count=evaluator.invalid_count, cache_hits=evaluator.cache_hits,
        elapsed_seconds=time.perf_counter() - t0, stop_reason=stop_reason,
        history=history,
    )


def greedy_hill_climb(
    budget: int = 300,
    seed: int = 42,
    max_proposal_factor: int = 100,
) -> RunResult:
    """贪心爬山：枚举当前序列的单点可变位点替换，采用 best-improvement。

    每轮将 FREE 位点及替代氨基酸顺序随机打乱（由 seed 控制），在预算内评估
    邻居；若本轮找到严格更优的合法邻居，则移动到该邻居继续搜索，否则停止。
    由于硬约束可能拒绝大量候选，实际 eval_count 可能低于预算；这是正常的。
    """
    if budget < 1:
        raise ValueError("budget 必须 >= 1")
    rng = random.Random(seed)
    evaluator = Evaluator()
    t0 = time.perf_counter()

    ref = evaluator.evaluate(REF_SEQ)
    assert ref is not None, "模块3参考序列未通过自身约束，请先修复模块3数据/配置"
    current = best = ref
    history = [_history_row("greedy_hill_climb", budget, seed, evaluator.eval_count,
                            evaluator.proposal_count, current, best, "reference")]
    visited = {REF_SEQ}
    max_proposals = max(budget * max_proposal_factor, budget + 10)
    stop_reason = "local_optimum"

    while evaluator.eval_count < budget and evaluator.proposal_count < max_proposals:
        positions = list(FREE)
        rng.shuffle(positions)
        best_neighbor: Optional[Evaluation] = None
        sweep_complete = True

        for pos in positions:
            alternatives = [aa for aa in ALPHABET if aa != current.sequence[pos]]
            rng.shuffle(alternatives)
            for aa in alternatives:
                if evaluator.eval_count >= budget or evaluator.proposal_count >= max_proposals:
                    sweep_complete = False
                    break
                chars = list(current.sequence)
                chars[pos] = aa
                cand = "".join(chars)
                if cand in visited:
                    continue
                visited.add(cand)
                ev = evaluator.evaluate(cand)
                if ev is None:
                    continue
                if ev.score > best.score:
                    best = ev
                if ev.score > current.score and (
                    best_neighbor is None or ev.score > best_neighbor.score
                ):
                    best_neighbor = ev
                history.append(_history_row(
                    "greedy_hill_climb", budget, seed, evaluator.eval_count,
                    evaluator.proposal_count, current, best))
            if not sweep_complete:
                break

        if best_neighbor is not None:
            current = best_neighbor
            history.append(_history_row(
                "greedy_hill_climb", budget, seed, evaluator.eval_count,
                evaluator.proposal_count, current, best, "move_to_better_neighbor"))
            if not sweep_complete:
                stop_reason = "evaluation_budget" if evaluator.eval_count >= budget else "proposal_limit"
                break
            continue

        if not sweep_complete:
            stop_reason = "evaluation_budget" if evaluator.eval_count >= budget else "proposal_limit"
        else:
            stop_reason = "local_optimum"
        break

    if evaluator.eval_count >= budget:
        stop_reason = "evaluation_budget"
    elif evaluator.proposal_count >= max_proposals:
        stop_reason = "proposal_limit"

    return RunResult(
        strategy="greedy_hill_climb", budget=budget, seed=seed, ref_id=REF_ID,
        reference_sequence=REF_SEQ, best_sequence=best.sequence,
        best_score=best.score, reference_score=ref.score,
        eval_count=evaluator.eval_count, proposal_count=evaluator.proposal_count,
        invalid_count=evaluator.invalid_count, cache_hits=evaluator.cache_hits,
        elapsed_seconds=time.perf_counter() - t0, stop_reason=stop_reason,
        history=history,
    )


def write_run(result: RunResult, outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "history.csv"), "w", newline="", encoding="utf-8-sig") as f:
        fields = list(result.history[0].keys()) if result.history else [
            "strategy", "budget", "seed", "eval", "proposal", "current_score",
            "best_so_far", "best_aggregation_max", "best_n_edits", "note"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(result.history)

    with open(os.path.join(outdir, "best.fasta"), "w", encoding="utf-8") as f:
        f.write(f">best|strategy={result.strategy}|budget={result.budget}|seed={result.seed}\n")
        f.write(result.best_sequence + "\n")
    with open(os.path.join(outdir, "reference.fasta"), "w", encoding="utf-8") as f:
        f.write(f">{result.ref_id}|reference\n{result.reference_sequence}\n")

    summary = asdict(result)
    summary.pop("history", None)
    with open(os.path.join(outdir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    config = {
        "strategy": result.strategy, "budget": result.budget, "seed": result.seed,
        "ref_id": result.ref_id, "reference_length": len(result.reference_sequence),
        "locked_sites_count": len(MASK), "free_sites_count": len(FREE),
        "objective": "maximize -aggregation_max; equivalent to minimize AGGRESCAN a3v 7-residue-window maximum",
        "hard_constraints": [
            "fixed sequence length and 20 standard amino acids",
            "locked positions unchanged",
            "natural property envelope from module3 profile (target aggregation_max excluded)",
        ],
        "pfam_validation": "PF00741 HMMER scan is performed after search on final best candidates, not in inner loop",
        "budget_definition": "maximum successful feasible objective evaluations including reference; infeasible proposals do not consume eval_count",
    }
    with open(os.path.join(outdir, "config.json"), "w", encoding="utf-8") as f:
        json.dump(config, f, ensure_ascii=False, indent=2)


def run_baseline(strategy: str, budget: int, seed: int, outdir: Optional[str] = None, **kwargs):
    if strategy == "random_mutation":
        result = random_mutation_search(budget=budget, seed=seed, **kwargs)
    elif strategy == "greedy_hill_climb":
        result = greedy_hill_climb(budget=budget, seed=seed, **kwargs)
    else:
        raise ValueError(f"未知策略: {strategy}")
    if outdir:
        write_run(result, outdir)
    return result
