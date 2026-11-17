"""成员3/B：GV02-04 定向优化的目标函数与约束

== 最终设计 ==
    max_s  target(s)        target = -聚集倾向(AGGRESCAN a3v 均值)
    s.t.   Pfam00741 通过 + 锁定位点不变 + 全部理化性质落在天然包络内

约束是**前提条件**，不是加权项；违反约束的候选直接判为不可行，不参与排序。
实现见 models/constraints.py，性质计算见 models/properties.py。

== 为什么是"聚集倾向最小化"==
任务书原文列出的可选目标包含「结构稳定性（AlphaFold pLDDT / 聚集倾向）」。
选它的理由是可控地验证过：野生型聚集倾向 0.3344，天然 856 条的最低值 0.2441，
有 27% 的真实下降空间，且需约 3 步突变可达 —— 既不平凡，也不靠冲出天然范围取胜。

== 为什么不能用别的（三次失败的设计，均有实测证据）==

[设计 1] 加权和代理：0.6*保守 + 0.3*组成 + 0.1*新颖
    corr(保守, 新颖) = -0.961 —— 两者本质是同一信息的正负两面，实际只有 1 个信号。
    更致命：最优解 = 家族共识序列。实测逐位最优把 conserv 推到 0.7063，
    高于 856 条天然序列的最好值 0.659 —— 即向平均值回归，不是设计。

[设计 2] 两亲性最大化：signed_amphiphilicity = 疏水性矩在野生型方向的投影
    爬山 15 步全部选 I（异亮氨酸），target +153%，但产出
    MAVEIVISISISSLAEVVDRILDKGIVIDIIVIISIVGIELIIIEARIVIASVITILIYAEAVGLISIAIVISV
    其性质全部冲出天然范围之外的极端值：
        GRAVY   0.775 -> 1.975（天然上限 0.867）
        聚集     0.334 -> 0.590（天然上限 0.363）
        pI      4.37  -> 3.88 （天然下限 3.99）
    幅度型指标（|μH|）方向无关，把强亲水残基放到同一相位也能抬高它；
    改投影虽修正了方向（最优突变转为疏水残基、Arg 被惩罚），却挡不住
    "整条序列都变疏水"这条退化路径。

[设计 3] 归一化两亲性（除以平均 |H| 使尺度不变）
    同样退化：聚集 0.535、GRAVY 1.683，仍远超天然上限。
    根因：天然序列本身已在各项理化性质的窄区间内，任何单向最大化
    都会冲出天然范围。

== 结论 ==
物理化学性质不能单独当目标，必须配「天然包络」硬约束（856 条天然序列实测的
取值范围，数据驱动、无主观阈值）。包络下 GRAVY 只剩 +0.09 余地，两亲性目标
基本无空间；而聚集倾向有 -0.09 余地，故取后者。

== 分工 ==
优化循环内调 score()            —— target + 约束，0.1 ms 级
候选产出后调 verify()           —— 追加 Pfam 官方模型扫描，秒级

== 面向优化器的接口（模块三交付）==
    GvpObjective            封装主目标的调用计数与可微性声明
        .is_differentiable  False（无梯度方法专用）
        .eval_count         成功通过约束的真实评估次数
        __call__(seq)       标量分值，违规抛 ConstraintViolation
        .feasible_score()   不可行返回 -inf，不计入 eval_count
        .evaluate(seq)      全部指标字典，含 eval_count
        .reset()            计数器归零
    default_objective       模块级单例；模块级 score/feasible_score/objective
                            均委托它计数，便于 RQ2 统计样本效率。
模块级函数签名与返回字段向后兼容（objective() 新增 eval_count 键）。
"""
import math
import os

try:
    from .constraints import (REF_ID, REF_SEQ, MASK, FREE, L, ConstraintViolation,
                              assert_ok, check, pfam_scan)
    from .properties import aggregation, all_properties, hydrophobic_moment
except ImportError:  # 直接以脚本运行
    from constraints import (REF_ID, REF_SEQ, MASK, FREE, L, ConstraintViolation,
                             assert_ok, check, pfam_scan)
    from properties import aggregation, all_properties, hydrophobic_moment

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

import json  # noqa: E402

DATA = os.environ.get(
    "GVP_DATA",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data"))
with open(os.path.join(DATA, "conservation_profile.json")) as _f:
    _P = json.load(_f)
_LOGP = [[math.log((_P["profile"][i][j] * _P["n_seq"] + 1.0) / (_P["n_seq"] + 20))
          for j in range(20)] for i in range(L)]
_ALIGNED = _P["aligned"]
ENVELOPE = _P["envelope"]


def target(seq):
    """主目标：越大越好，等价于聚集倾向越小。

    取 aggregation_max（最长 7 残基窗口的 a3v 均值）而非 aggregation_mean：
    mean 是逐残基可加的 Σa(s_i)/L，对残基位置不敏感，其无约束最小值只在
    组成上取到（全脯氨酸，实测 agg_mean=-0.4696 越界）—— 优化问题退化成
    选组成。max 依赖残基位置，优化器必须定位最强的聚集热点，搜索才非平凡。
    """
    return -aggregation(seq)[1]


def conserv(seq):
    """辅助指标：位点保守性（各残基天然概率的几何平均）。诊断用，不进目标。"""
    s = 0.0
    for i, aa in enumerate(seq):
        s += _LOGP[i][ALPHABET.index(aa)]
    return math.exp(s / L)


def novelty(seq):
    """辅助指标：1 - 与天然集合最大一致度。诊断用，不进目标。"""
    best = max(sum(1 for i, aa in enumerate(seq) if row[i] == aa)
               for row in _ALIGNED)
    return 1.0 - best / L


def _score_raw(seq):
    """主目标原值。约束违规抛 ConstraintViolation。无计数副作用。"""
    assert_ok(seq)
    return target(seq)


def _objective_raw(seq):
    """主目标 + 全部辅助指标。违规抛 ConstraintViolation。无计数副作用。"""
    assert_ok(seq)
    out = all_properties(seq)
    out.update({
        "target": round(target(seq), 4),
        "aggregation_mean": round(aggregation(seq)[0], 4),
        "aggregation_max": round(aggregation(seq)[1], 4),
        "conserv": round(conserv(seq), 4),
        "novelty": round(novelty(seq), 4),
        "muH_beta": round(hydrophobic_moment(seq, 180.0), 4),
        "n_edits": sum(1 for a, b in zip(seq, REF_SEQ) if a != b),
    })
    return out


class GvpObjective:
    """面向优化器的目标函数封装。

    is_differentiable = False：目标是 max() 的负值、约束含离散判定，
    梯度不存在 → 优化器只能用无梯度方法（GA / BO）。

    eval_count 只统计**成功通过约束的真实评估**；不可行候选不计入，
    以保证 RQ2 的「样本效率」口径一致。
    """

    is_differentiable = False

    def __init__(self):
        self.eval_count = 0

    def __call__(self, seq):
        v = _score_raw(seq)
        self.eval_count += 1
        return v

    def evaluate(self, seq):
        out = _objective_raw(seq)
        self.eval_count += 1
        out["eval_count"] = self.eval_count
        return out

    def feasible_score(self, seq, infeasible=float("-inf")):
        """不可行时返回 infeasible 且不计入 eval_count，不抛异常。

        只捕获 ConstraintViolation —— 这是唯一预期的不可行信号。
        其余异常（KeyError / TypeError / 其他 ValueError）属内部缺陷，
        必须向上冒泡，不得伪装成「不可行候选」。
        """
        try:
            v = _score_raw(seq)
        except ConstraintViolation:
            return infeasible
        self.eval_count += 1
        return v

    def reset(self):
        self.eval_count = 0


default_objective = GvpObjective()


def score(seq):
    """优化循环入口，计入 default_objective.eval_count。违规抛 ConstraintViolation。"""
    return default_objective(seq)


def feasible_score(seq, infeasible=float("-inf")):
    """优化器安全包装：不可行返回 infeasible 且不计入评估数，不抛异常。"""
    return default_objective.feasible_score(seq, infeasible)


def objective(seq):
    """完整评估（含 eval_count 键），计入 default_objective.eval_count。"""
    return default_objective.evaluate(seq)


def verify(seqs, names=None):
    """候选产出后的重验证：Pfam00741 官方模型 + gathering threshold。"""
    seqs = list(seqs)
    names = list(names) if names else [f"cand{i}" for i in range(len(seqs))]
    return pfam_scan(seqs, names)


def greedy(seq=None, max_steps=15, verbose=False):
    """简单爬山基线（M3 用）。每步取提升最大且仍可行的单点突变。"""
    cur = seq or REF_SEQ
    assert_ok(cur)
    for step in range(max_steps):
        b = target(cur)
        bp, ba, bv = None, None, b
        for p in FREE:
            for a in ALPHABET:
                if a == cur[p]:
                    continue
                v = feasible_score(cur[:p] + a + cur[p + 1:])
                if v > bv:
                    bp, ba, bv = p, a, v
        if bp is None:
            break
        cur = cur[:bp] + ba + cur[bp + 1:]
        if verbose:
            print(f"  step {step+1}: {bp} {REF_SEQ[bp]}->{ba}  {b:.4f} -> {bv:.4f}")
    return cur


def _selftest():
    import random
    import time

    w = objective(REF_SEQ)
    w2 = objective(REF_SEQ)
    assert ({k: v for k, v in w.items() if k != "eval_count"} ==
            {k: v for k, v in w2.items() if k != "eval_count"}), "非确定性"
    print(f"[1] 决定性 OK  参考 {REF_ID} len={L}")
    print(f"    野生型 聚集热点={w['aggregation_max']:.4f} (均值={w['aggregation_mean']:.4f}) "
          f"conserv={w['conserv']:.4f} GRAVY={w['gravy']:.4f} pI={w['pI']}")
    print(f"    天然包络 聚集均值=[{ENVELOPE['aggregation_mean']['min']:.4f},"
          f"{ENVELOPE['aggregation_mean']['max']:.4f}]  "
          f"GRAVY=[{ENVELOPE['gravy']['min']:.4f},{ENVELOPE['gravy']['max']:.4f}]")
    print(f"    目标性质 aggregation_max 不入包络约束（否则极限恒等于约束边界）")

    for bad, why in ((REF_SEQ + "A", "长度"), (REF_SEQ[:5] + "X" + REF_SEQ[6:], "非法残基"),
                     (REF_SEQ[:MASK[0]] + "W" + REF_SEQ[MASK[0] + 1:], "锁定改动")):
        try:
            score(bad)
        except ConstraintViolation:
            continue
        raise AssertionError(f"{why} 违规未拦截")
    print("[2] 约束违规 fail-fast OK")

    # 包络约束必须能**单独**拦住"结构合法但理化性质冲出天然范围"的序列。
    # 用 FREE 位点全替换为 Ile 构造：MASK 不动 → 结构合法，只违反包络。
    # 不能用设计 2 的爬山输出：那条序列同时改了 32 个锁定位点，是被结构约束
    # 拦下的，无法证明包络在起作用。
    _free = set(FREE)
    env_bad = "".join("I" if i in _free else REF_SEQ[i] for i in range(L))
    _v = check(env_bad)
    assert _v and not any(("锁定" in x or "长度" in x or "非法" in x) for x in _v), \
        f"该序列应只违反包络，实际 {_v}"
    assert feasible_score(env_bad) == float("-inf"), "天然包络未拦住越界序列"
    print(f"[3] 天然包络单独拦住性质越界序列 OK（结构合法，{len(_v)} 项越界）")

    random.seed(0)
    feas, allm = [], []
    for p in FREE:
        for a in ALPHABET:
            if a == REF_SEQ[p]:
                continue
            s = REF_SEQ[:p] + a + REF_SEQ[p + 1:]
            allm.append(target(s))
            if check(s) == []:
                feas.append(target(s))
    print(f"[4] 单点突变 {len(allm)} 条，其中可行 {len(feas)} 条"
          f"（{len(feas)/len(allm)*100:.0f}%）")
    if feas:
        print(f"    可行区间 target ∈ [{min(feas):.4f}, {max(feas):.4f}]"
              f" 跨度 {max(feas)-min(feas):.4f}")

    t = time.time()
    for _ in range(200):
        feasible_score(REF_SEQ)
    print(f"[5] score+约束 单次 {(time.time()-t)*5:.3f} ms")

    print("[6] 爬山基线（受限于约束）...")
    g = greedy(verbose=True)
    o = objective(g)
    print(f"    改动 {o['n_edits']} 位")
    print(f"    聚集热点 {w['aggregation_max']:.4f} -> {o['aggregation_max']:.4f}"
          f"  ({(o['aggregation_max']-w['aggregation_max'])/w['aggregation_max']*100:+.1f}%)"
          f"   天然最低 {ENVELOPE['aggregation_max']['min']:.4f}")
    print(f"    聚集均值 {w['aggregation_mean']:.4f} -> {o['aggregation_mean']:.4f}")
    print(f"    conserv {w['conserv']:.4f} -> {o['conserv']:.4f}   "
          f"GRAVY {w['gravy']:.4f} -> {o['gravy']:.4f}   pI {w['pI']} -> {o['pI']}")
    from constraints import ENVELOPE_SKIP
    checked = [k for k in ENVELOPE if k not in ENVELOPE_SKIP]
    for k in checked:
        val = o[k]
        assert ENVELOPE[k]["min"] <= val <= ENVELOPE[k]["max"], f"{k} 越界"
    print(f"    受约束的 {len(checked)} 项理化性质仍在天然包络内 OK"
          f"（目标项 {list(ENVELOPE_SKIP)} 除外）")
    print(f"    距离各项包络边界的余量: " + "  ".join(
        f"{k}={min(o[k]-ENVELOPE[k]['min'], ENVELOPE[k]['max']-o[k]):.4f}" for k in checked))

    r = verify([REF_SEQ, g], ["wildtype", "greedy"])
    for x in r:
        print(f"    {x['name']:9} Pfam00741 E={x['evalue']:.2e} passed={x['passed']}")
    assert all(x["passed"] for x in r), "优化后未通过 Pfam"
    print("[7] 优化后仍通过 Pfam00741 官方模型 OK")

    # ---- GvpObjective 封装自检 ----
    obj = GvpObjective()
    assert GvpObjective.is_differentiable is False and obj.is_differentiable is False
    assert obj.eval_count == 0
    v0 = obj(REF_SEQ)
    assert isinstance(v0, float) and abs(v0 - target(REF_SEQ)) < 1e-12
    assert obj.eval_count == 1, "计次未生效"
    for _ in range(4):
        obj(REF_SEQ)
    assert obj.eval_count == 5, f"计次错误 {obj.eval_count}"
    d = obj.evaluate(REF_SEQ)
    need = {"target", "aggregation_max", "aggregation_mean", "muH_beta", "muH_helix",
            "beta_propensity", "gravy", "pI", "net_charge_pH7", "conserv",
            "novelty", "n_edits", "eval_count"}
    assert need <= set(d), f"evaluate 缺字段 {need - set(d)}"
    assert d["eval_count"] == 6
    assert obj.feasible_score(env_bad) == float("-inf")
    assert obj.eval_count == 6, "不可行候选不应计入评估数"
    obj.reset()
    assert obj.eval_count == 0
    print(f"[8] GvpObjective OK  is_differentiable=False  eval_count 计次/重置/不可行不计次")
    print("\n全部自检通过")


if __name__ == "__main__":
    _selftest()
