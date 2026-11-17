"""硬约束检查（GV02-04）

任务书原文：「在保持 Pfam 域约束满足和关键位点保守的前提下，定向优化某项目标
性质」—— 约束是**前提条件**，不是目标函数里的加权项。因此这里独立成层：

    Tier 1（纯 python，微秒级）→ 优化循环内每次调用
        check(seq)       全部硬约束，返回违规描述列表
        satisfies(seq)   布尔可行性快判
        penalty(seq)     违规程度的连续惩罚量（可行序列为 0.0）
        project(seq)     投影回结构可行子空间 C_struct
    Tier 2（HMMER 子进程，秒级）→ **只对最终候选集批量跑，禁止进循环**
        pfam_scan(seqs)  官方 Pfam 模型扫描

快速约束（违反即抛异常，不静默通过）：
    1. 长度必须等于参考长度（本优化为 substitution-only）
    2. 只能是 20 种标准氨基酸
    3. MASK 锁定位点必须与参考序列一致
    4. 6 项理化性质落在 856 条天然序列实测包络内（目标项 aggregation_max 除外）

Pfam 约束用 PF00741（Gas vesicle protein）官方模型 + Pfam gathering threshold
(--cut_ga)。注意：任务书写的是 PF01132，但那是 Elongation factor P 的 OB 域，
扫全部 856 条天然 GvpA 零命中；PF00741 才对参考序列给出 E=1.5e-21 的命中。
详见 README。

面向优化器的封装见 GvpConstraints 与模块底单例 default_constraints。
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("GVP_DATA", os.path.join(HERE, "..", "..", "data"))
PFAM_DIR = os.path.join(DATA, "pfam")
PFAM_HMM = os.environ.get("GVP_PFAM_HMM", os.path.join(PFAM_DIR, "PF00741.hmm"))

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
_AA_SET = frozenset(ALPHABET)

import json  # noqa: E402

_PROFILE_PATH = os.path.join(DATA, "conservation_profile.json")
if not os.path.exists(_PROFILE_PATH):
    raise RuntimeError(
        f"缺少保守性图谱 {_PROFILE_PATH}\n"
        f"请先运行: python3 preprocessing/build_conservation_profile.py")

with open(_PROFILE_PATH) as _f:
    _P = json.load(_f)

REF_ID = _P["ref_id"]
REF_SEQ = _P["ref_seq"]
L = len(REF_SEQ)
MASK = tuple(_P["mask"])
FREE = tuple(i for i in range(L) if i not in set(MASK))


class ConstraintViolation(ValueError):
    """违反硬约束。优化循环应把它当作「该候选不可行」，而非静默通过。"""


# 天然性质包络：856 条天然 GvpA 实测的取值范围。
# 为什么需要它：实测表明单向最大化任何理化性质都会冲出天然范围 —— 例如
# 爬山基线把 target 从 0.0236 抬到 0.0598，代价是 GRAVY 0.775->1.975
#（天然上限 0.867）、聚集 0.334->0.590（天然上限 0.363），产出一条高疏水、
# 易聚集的非天然序列。所有阈值均取自数据，不设主观数值。
_ENV = _P.get("envelope", {})
ENVELOPE_KEYS = tuple(_ENV.keys())
# "minmax" (默认，最宽) 或 "p5p95"（排除尾部，更严）
ENVELOPE_MODE = os.environ.get("GVP_ENVELOPE_MODE", "minmax")

# 目标性质本身不进包络约束 —— 既约束又优化同一个量会使最优解恒等于约束边界，
# 使 RQ1「提升极限」退化成一个恒等式。剔除后，极限由其余约束的交织决定。
ENVELOPE_SKIP = tuple(
    s.strip() for s in os.environ.get("GVP_ENVELOPE_SKIP", "aggregation_max").split(",")
    if s.strip())


def check(seq):
    """全部硬约束。返回空列表表示通过，否则返回违规描述列表。"""
    v = []
    if not isinstance(seq, str):
        return [f"非字符串输入: {type(seq).__name__}"]
    if len(seq) != L:
        v.append(f"长度 {len(seq)} != 参考长度 {L}（本优化为 substitution-only）")
    bad = set(seq) - set(ALPHABET)
    if bad:
        v.append(f"非法残基 {sorted(bad)}")
    if len(seq) == L and not bad:
        diff = [i for i in MASK if seq[i] != REF_SEQ[i]]
        if diff:
            v.append(f"锁定位点被改动 {len(diff)} 处: {diff[:10]}")
    if v:
        return v          # 结构已非法，性质计算无意义
    return check_properties(seq)


def assert_ok(seq):
    """通过则返回 None，否则抛 ConstraintViolation。"""
    v = check(seq)
    if v:
        raise ConstraintViolation("; ".join(v))


def _property_values(seq):
    """当前序列的 6 项包络属性实测值。check_properties 与 penalty 共用。"""
    try:
        from .properties import (aggregation, beta_propensity, gravy as _g,
                                 isoelectric_point, net_charge)
    except ImportError:
        from properties import (aggregation, beta_propensity, gravy as _g,
                                isoelectric_point, net_charge)
    am, ax = aggregation(seq)
    return {"aggregation_mean": am, "aggregation_max": ax, "gravy": _g(seq),
            "beta_propensity": beta_propensity(seq),
            "pI": isoelectric_point(seq), "net_charge_pH7": net_charge(seq, 7.0)}


def _envelope_bounds(k):
    """给定属性的包络上下界；不在包络内或被剔除则返回 None。"""
    env = _ENV.get(k)
    if not env or k in ENVELOPE_SKIP:
        return None
    lo_key, hi_key = ("min", "max") if ENVELOPE_MODE == "minmax" else ("p5", "p95")
    return env[lo_key], env[hi_key]


def check_properties(seq):
    """天然性质包络检查。返回违规描述列表。"""
    v = []
    for k, val in _property_values(seq).items():
        b = _envelope_bounds(k)
        if b is None:
            continue
        lo, hi = b
        if not (lo <= val <= hi):
            v.append(f"{k}={val:.4f} 超出天然范围 [{lo:.4f}, {hi:.4f}]")
    return v


def satisfies(seq):
    """Tier 1 快速可行性判定，等价于 check(seq) 为空。纯 python，微秒级。"""
    return not check(seq)


def project(seq):
    """投影到结构可行子空间 C_struct：长度 L、合法字母表、MASK 位点复原。

    仅处理结构约束。天然包络 C_env 是 5 维非凸区间，无确定性单步投影，
    其违规量由 penalty() 给出。
    输出必定满足 C_struct（长度 L、字母全合法、锁定位点 100% 等于参考）。
    """
    if not isinstance(seq, str):
        raise TypeError(f"序列必须是 str，收到 {type(seq).__name__}")
    chars = list(seq[:L])
    if len(chars) < L:
        chars.extend(REF_SEQ[len(chars):L])
    for i, a in enumerate(chars):
        if a not in _AA_SET:
            chars[i] = REF_SEQ[i]
    for i in MASK:
        chars[i] = REF_SEQ[i]
    return "".join(chars)


# 惩罚权重。仅供软约束/罚函数类优化器；硬约束场景用 satisfies。
PEN_LEN, PEN_ILLEGAL, PEN_MASK, PEN_ENVELOPE = 10.0, 10.0, 5.0, 2.0


def penalty(seq):
    """违反约束的连续惩罚量，可行序列严格返回 0.0。

    结构违规与包络越界不叠加：结构非法时性质计算无意义，直接返回结构惩罚。
    包络项按区间宽度归一化为违规比例。
    """
    if not isinstance(seq, str):
        return PEN_LEN * (L + 1)
    p = PEN_LEN * abs(len(seq) - L)
    p += PEN_ILLEGAL * sum(1 for a in seq if a not in _AA_SET)
    n = min(len(seq), L)
    p += PEN_MASK * sum(1 for i in MASK if i < n and seq[i] != REF_SEQ[i])
    if p > 0.0:
        return p
    for k, val in _property_values(seq).items():
        b = _envelope_bounds(k)
        if b is None:
            continue
        lo, hi = b
        w = hi - lo
        if w <= 0.0:
            continue
        if val < lo:
            p += PEN_ENVELOPE * (lo - val) / w
        elif val > hi:
            p += PEN_ENVELOPE * (val - hi) / w
    return p


# ---------------------------------------------------------------- Pfam 扫描

def _hmmer():
    exe = os.environ.get("HMMSEARCH_BIN") or shutil.which("hmmsearch")
    if not exe:
        raise RuntimeError(
            "找不到 hmmsearch。安装：sudo apt-get install hmmer\n"
            "或设置 HMMSEARCH_BIN 环境变量")
    return exe


def pfam_scan(seqs, names=None, hmm=None):
    """用官方 Pfam 模型 + gathering threshold 扫描。

    seqs: 序列字符串可迭代；names: 对应名称
    返回: [{"name","evalue","score","passed"}]，passed 为 --cut_ga 判定结果
    """
    hmm = hmm or PFAM_HMM
    if not os.path.exists(hmm):
        raise RuntimeError(f"找不到 Pfam 模型 {hmm}")
    seqs = list(seqs)
    names = list(names) if names else [f"seq{i}" for i in range(len(seqs))]
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "in.fa")
        out = os.path.join(td, "out.domtbl")
        with open(fa, "w") as f:
            for n, s in zip(names, seqs):
                f.write(f">{n}\n{s}\n")
        subprocess.run(
            [_hmmer(), "--noali", "--cut_ga", "--domtblout", out, hmm, fa],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        hits = {}
        with open(out) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                p = line.split()
                if len(p) < 8:
                    continue
                # 取该序列的最佳 domain（domtbl 按 score 降序，首次出现即最佳）
                hits.setdefault(p[0], (float(p[6]), float(p[7])))
    return [{"name": n,
             "evalue": hits.get(n, (float("inf"), 0.0))[0],
             "score": hits.get(n, (float("inf"), 0.0))[1],
             "passed": n in hits}
            for n in names]


class GvpConstraints:
    """面向优化器的约束层封装。

    Tier 1（satisfies/project/penalty/check/assert_ok）为纯 python，微秒级，
    可在优化循环内每次调用。
    Tier 2（verify）会起 HMMER 子进程，秒级，**只能对最终候选集批量调用**。

    这些方法是对模块级函数的稳定门面（facade），不是冗余抽象：数据集只有
    唯一一份（由 conservation_profile.json 加载），因此不可能存在持有不同
    数据的第二个实例，门面只为让模块四拿到一个可传递的句柄与一组固定入口。
    """

    def __init__(self):
        self.ref_id = REF_ID
        self.ref_seq = REF_SEQ
        self.mask = MASK
        self.free = FREE
        self.envelope = _ENV
        self.envelope_skip = ENVELOPE_SKIP
        self.envelope_mode = ENVELOPE_MODE

    def satisfies(self, seq):
        return satisfies(seq)

    def project(self, seq):
        return project(seq)

    def penalty(self, seq):
        return penalty(seq)

    def check(self, seq):
        return check(seq)

    def assert_ok(self, seq):
        return assert_ok(seq)

    def verify(self, seqs, names=None):
        """Tier 2：PF00741 官方模型批量重判。禁止在优化循环内调用。"""
        return pfam_scan(seqs, names)


default_constraints = GvpConstraints()


def _selftest():
    print(f"参考 {REF_ID}  len={L}  锁定={len(MASK)}  可改={len(FREE)}")

    assert check(REF_SEQ) == [], "参考序列不该违规"
    assert check(REF_SEQ + "A"), "长度违规未检出"
    assert check(REF_SEQ[:5] + "X" + REF_SEQ[6:]), "非法残基未检出"
    p0 = MASK[0]
    assert check(REF_SEQ[:p0] + "A" + REF_SEQ[p0 + 1:]) or REF_SEQ[p0] == "A", \
        "锁定位点改动未检出"
    print("[1] 快速约束 OK（非字符串 / 长度 / 非法残基 / 锁定改动 全部拦截）")

    try:
        assert_ok(REF_SEQ[:p0] + "W" + REF_SEQ[p0 + 1:])
    except ConstraintViolation:
        pass
    else:
        raise AssertionError("assert_ok 未抛异常")
    print("[2] ConstraintViolation 抛错 OK")

    # satisfies / project / penalty
    assert satisfies(REF_SEQ) is True, "参考序列应满足约束"
    assert satisfies(REF_SEQ[:5] + "X" + REF_SEQ[6:]) is False, "非法残基未被判不可行"
    assert penalty(REF_SEQ) == 0.0, "可行序列惩罚必须为 0"
    assert penalty(REF_SEQ[:5] + "X" + REF_SEQ[6:]) > 0.0, "非法残基惩罚必须为正"
    pens = []
    for k in (1, 3, 5):
        s = list(REF_SEQ)
        for i in MASK[:k]:
            s[i] = "W" if REF_SEQ[i] != "W" else "A"
        pens.append(penalty("".join(s)))
    assert pens[0] < pens[1] < pens[2], f"惩罚未随违规数单调递增: {pens}"

    proj = project(REF_SEQ[:10] + "X" * 5 + REF_SEQ[15:])
    assert len(proj) == L and set(proj) <= _AA_SET, "project 未修复长度/字符"
    assert all(proj[i] == REF_SEQ[i] for i in MASK), "project 未复原锁定位点"
    assert len(project("M" * 30)) == L, "project 未补齐短序列"
    assert len(project("M" * 200)) == L, "project 未截断长序列"
    print(f"[3] satisfies/project/penalty OK  锁定位点惩罚 {pens}")

    c = default_constraints
    assert c.ref_seq == REF_SEQ and len(c.mask) == len(MASK) and len(c.free) == len(FREE)
    assert c.satisfies(REF_SEQ) and c.penalty(REF_SEQ) == 0.0
    assert c.project("M" * 5) == project("M" * 5)
    print(f"[4] GvpConstraints 封装 OK  参考={c.ref_id} 包络项={len(c.envelope)}")

    if shutil.which("hmmsearch") and os.path.exists(PFAM_HMM):
        r = pfam_scan([REF_SEQ, REF_SEQ[:5] + "W" * (L - 5)],
                      names=["ref", "garbage"])
        print(f"[5] Pfam 扫描  ref: E={r[0]['evalue']:.2e} passed={r[0]['passed']}")
        print(f"    garbage: passed={r[1]['passed']}")
        assert r[0]["passed"], "参考序列应通过 Pfam00741"
        assert not r[1]["passed"], "乱序列不应通过"
        assert c.verify([REF_SEQ], ["wt"])[0]["passed"], "GvpConstraints.verify 未通过"
        print("[6] 可通过参考、拦得住垃圾序列 / verify 封装 OK")
    else:
        print("[5] 跳过 Pfam 扫描（未安装 hmmer 或缺模型）")
    print("\n全部自检通过")


if __name__ == "__main__":
    _selftest()
