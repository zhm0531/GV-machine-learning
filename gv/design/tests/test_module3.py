"""模块三单元测试：目标函数与约束层接口契约。

运行（在 gv/design 下）：
    python3 -m unittest tests/test_module3.py -v
    python3 tests/test_module3.py

纯标准库，无第三方依赖。Tier 2（Pfam HMMER）在缺少 hmmsearch 或模型文件时跳过。
"""
import os
import shutil
import sys
import time
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, "..", "src"))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from models.constraints import (  # noqa: E402
    ALPHABET, ENVELOPE_SKIP, FREE, L, MASK, REF_SEQ,
    ConstraintViolation, GvpConstraints, default_constraints, pfam_scan,
)
from models.objective import ENVELOPE, GvpObjective, default_objective  # noqa: E402

# 只违反包络、结构完全合法的序列：FREE 位点全替换为 Ile（MASK 不动）。
# 不能用设计 2 的爬山输出：那条序列同时改了 32 个锁定位点，会被结构约束拦下，
# 无法证明包络约束在起作用。
_FREE = set(FREE)
ENV_BAD = "".join("I" if i in _FREE else REF_SEQ[i] for i in range(L))

ILLEGAL = REF_SEQ[:5] + "X" + REF_SEQ[6:]
TOO_LONG = REF_SEQ + "A"
MASK_HIT = REF_SEQ[:MASK[0]] + ("W" if REF_SEQ[MASK[0]] != "W" else "A") \
    + REF_SEQ[MASK[0] + 1:]


class TestGvpConstraints(unittest.TestCase):
    """Tier 1 约束接口 + Tier 2 Pfam 重判。"""

    def setUp(self):
        self.c = GvpConstraints()

    def test_reference_satisfies(self):
        self.assertTrue(self.c.satisfies(REF_SEQ))
        self.assertEqual(self.c.check(REF_SEQ), [])
        self.assertEqual(self.c.penalty(REF_SEQ), 0.0)

    def test_structural_violations(self):
        for bad, why in ((TOO_LONG, "长度"), (ILLEGAL, "非法残基"), (MASK_HIT, "锁定改动")):
            with self.subTest(why=why):
                self.assertFalse(self.c.satisfies(bad))
                self.assertTrue(self.c.check(bad))
                self.assertGreater(self.c.penalty(bad), 0.0)
                with self.assertRaises(ConstraintViolation):
                    self.c.assert_ok(bad)

    def test_envelope_violations(self):
        # 结构合法（无长度/字符/掩码违规），仅包络越界
        v = self.c.check(ENV_BAD)
        self.assertTrue(v, "该序列应违反包络")
        self.assertFalse(any(("长度" in x or "非法" in x or "锁定" in x) for x in v), v)
        self.assertTrue(any("超出天然范围" in x for x in v), v)
        self.assertFalse(self.c.satisfies(ENV_BAD))
        self.assertGreater(self.c.penalty(ENV_BAD), 0.0)
        # 目标性质本身不入包络
        self.assertIn("aggregation_max", ENVELOPE_SKIP)
        self.assertFalse(any("aggregation_max" in x for x in v), v)

    def test_project_operator(self):
        victim = list(REF_SEQ)
        for i in list(MASK[:5]) + list(FREE[:5]):
            victim[i] = "X" if i in MASK else "W"
        p = self.c.project("".join(victim))
        self.assertEqual(len(p), L)
        self.assertTrue(set(p) <= set(ALPHABET))
        for i in MASK:
            self.assertEqual(p[i], REF_SEQ[i])
        # 长度异常与非法字符一律修复
        self.assertEqual(len(self.c.project("M" * 5)), L)
        self.assertEqual(len(self.c.project("M" * 300)), L)
        self.assertTrue(set(self.c.project("Z" * L)) <= set(ALPHABET))
        with self.assertRaises(TypeError):
            self.c.project(None)

    def test_penalty_monotonicity(self):
        pens = []
        for k in (1, 3, 5):
            s = list(REF_SEQ)
            for i in MASK[:k]:
                s[i] = "W" if REF_SEQ[i] != "W" else "A"
            pens.append(self.c.penalty("".join(s)))
        self.assertEqual(pens, [5.0, 15.0, 25.0])
        self.assertTrue(pens[0] < pens[1] < pens[2])

    @unittest.skipUnless(shutil.which("hmmsearch"), "未安装 hmmsearch，跳过 Tier 2")
    def test_pfam_two_tier(self):
        r = self.c.verify([REF_SEQ, ENV_BAD], ["wt", "env_bad"])
        self.assertEqual([x["name"] for x in r], ["wt", "env_bad"])
        self.assertTrue(r[0]["passed"], f"参考序列应通过 PF00741: {r[0]}")
        # E-value 随搜索空间大小变化（同序列单条搜索 1.2e-21，两条 2.4e-21），
        # 故只比对与批次无关的 score / passed。
        single = pfam_scan([REF_SEQ], ["wt"])[0]
        self.assertEqual(r[0]["score"], single["score"])
        self.assertEqual(r[0]["passed"], single["passed"])
        self.assertEqual(set(r[0]), {"name", "evalue", "score", "passed"})
        # 两级开销必须分离：Tier 1 纯 python，Tier 2 要起 HMMER 子进程
        t = time.perf_counter()
        for _ in range(20):
            self.c.satisfies(REF_SEQ)
        t_sat = (time.perf_counter() - t) / 20
        t = time.perf_counter()
        pfam_scan([REF_SEQ], ["wt"])
        t_hmm = time.perf_counter() - t
        self.assertLess(t_sat * 5, t_hmm,
                        f"Tier1 应远快于 Tier2：{t_sat*1000:.3f} ms vs {t_hmm*1000:.2f} ms")

    def test_tier1_is_fast(self):
        """satisfies 是毫秒以下快判，可在优化循环内每次调用（不依赖 HMMER）。"""
        t = time.perf_counter()
        for _ in range(200):
            self.c.satisfies(REF_SEQ)
        per_call = (time.perf_counter() - t) / 200
        self.assertLess(per_call, 2e-3, f"实测 {per_call*1000:.3f} ms/次")


class TestGvpObjective(unittest.TestCase):
    """目标函数接口契约。"""

    def setUp(self):
        self.obj = GvpObjective()

    def test_callable_interface(self):
        w = default_objective.evaluate(REF_SEQ)
        self.assertAlmostEqual(w["aggregation_max"], 0.7134, places=3)
        v = self.obj(REF_SEQ)
        self.assertIsInstance(v, float)
        # 规格写明的期望值：-φ_max(野生型) = -0.7134（四位小数）
        self.assertAlmostEqual(round(v, 4), -0.7134, places=6)
        # target() 返回未舍入值，evaluate() 的字段舍入到 4 位
        self.assertAlmostEqual(v, -w["aggregation_max"], places=4)

    def test_is_differentiable(self):
        self.assertIs(GvpObjective.is_differentiable, False)
        self.assertIs(self.obj.is_differentiable, False)
        self.assertIs(default_objective.is_differentiable, False)

    def test_eval_count_tracking(self):
        self.assertEqual(self.obj.eval_count, 0)
        for _ in range(5):
            self.obj(REF_SEQ)
        self.assertEqual(self.obj.eval_count, 5)
        self.obj.reset()
        self.assertEqual(self.obj.eval_count, 0)

    def test_evaluate_schema(self):
        d = self.obj.evaluate(REF_SEQ)
        need = {"target", "aggregation_max", "aggregation_mean", "muH_beta", "muH_helix",
                "beta_propensity", "gravy", "pI", "net_charge_pH7", "conserv",
                "novelty", "n_edits", "eval_count"}
        self.assertTrue(need <= set(d), f"缺字段 {need - set(d)}")
        self.assertEqual(d["n_edits"], 0)
        self.assertEqual(d["eval_count"], 1)
        self.assertAlmostEqual(d["target"], -d["aggregation_max"], places=6)
        for k in ("target", "aggregation_max", "aggregation_mean", "muH_beta",
                  "muH_helix", "beta_propensity", "gravy", "pI",
                  "net_charge_pH7", "conserv", "novelty"):
            self.assertIsInstance(d[k], float, k)
        for k in ("n_edits", "eval_count"):
            self.assertIsInstance(d[k], int, k)
        # 取值范围
        self.assertGreaterEqual(d["conserv"], 0.0)
        self.assertLessEqual(d["conserv"], 1.0)
        self.assertGreaterEqual(d["novelty"], 0.0)
        self.assertLessEqual(d["novelty"], 1.0)
        self.assertGreaterEqual(d["pI"], 0.0)
        self.assertLessEqual(d["pI"], 14.0)
        self.assertGreaterEqual(d["muH_beta"], 0.0)
        self.assertGreaterEqual(d["muH_helix"], 0.0)
        self.assertGreaterEqual(d["beta_propensity"], 0.0)
        # a3v 量表取值范围 [-1.115, 1.161]，留宽松余量做区间哨兵
        for k in ("aggregation_max", "aggregation_mean"):
            self.assertGreaterEqual(d[k], -2.0, k)
            self.assertLessEqual(d[k], 2.0, k)
        env = ENVELOPE
        for k in env:
            if k in ENVELOPE_SKIP:
                continue
            self.assertLessEqual(env[k]["min"], d[k])
            self.assertLessEqual(d[k], env[k]["max"])

    def test_internal_error_propagates(self):
        """非 ConstraintViolation 的内部缺陷必须冒泡，不得伪装成「不可行候选」。"""
        import models.objective as O
        with mock.patch.object(O, "_score_raw", side_effect=KeyError("boom")):
            with self.assertRaises(KeyError):
                self.obj.feasible_score(REF_SEQ)
        self.assertEqual(self.obj.eval_count, 0)

    def test_constraint_fail_fast(self):
        for bad in (TOO_LONG, ILLEGAL, MASK_HIT, ENV_BAD):
            with self.subTest(bad=bad[:12]):
                with self.assertRaises(ConstraintViolation):
                    self.obj(bad)
                with self.assertRaises(ConstraintViolation):
                    self.obj.evaluate(bad)

    def test_feasible_score(self):
        self.assertEqual(self.obj.feasible_score(ENV_BAD), float("-inf"))
        self.assertEqual(self.obj.feasible_score(ILLEGAL), float("-inf"))
        self.assertEqual(self.obj.eval_count, 0, "不可行候选不应计入评估数")
        v = self.obj.feasible_score(REF_SEQ)
        self.assertGreater(v, float("-inf"))
        self.assertEqual(self.obj.eval_count, 1)
        self.assertEqual(self.obj.feasible_score(ENV_BAD, infeasible=-1e9), -1e9)

    def test_determinism(self):
        a = self.obj.evaluate(REF_SEQ)
        b = self.obj.evaluate(REF_SEQ)
        a.pop("eval_count"), b.pop("eval_count")
        self.assertEqual(a, b)
        self.assertEqual(self.obj(REF_SEQ), self.obj(REF_SEQ))


class TestModule3Integration(unittest.TestCase):
    """模拟模块四优化器的一轮步进循环。"""

    def test_optimizer_step_loop(self):
        import random

        cons = default_constraints
        obj = GvpObjective()
        random.seed(42)

        # 采样 20 条可行候选：单点突变，仅落在 FREE 位点
        cands = []
        while len(cands) < 20:
            p = random.choice(FREE)
            a = random.choice([x for x in ALPHABET if x != REF_SEQ[p]])
            c = cons.project(REF_SEQ[:p] + a + REF_SEQ[p + 1:])
            if cons.satisfies(c):
                cands.append(c)
        self.assertEqual(len(set(cands)), len(cands), "候选应互不相同")
        for c in cands:
            self.assertEqual(len(c), L)
            for i in MASK:
                self.assertEqual(c[i], REF_SEQ[i])
            self.assertEqual(cons.penalty(c), 0.0)

        self.assertEqual(obj.eval_count, 0)
        scored = [(obj(c), c) for c in cands]
        self.assertEqual(obj.eval_count, 20, "评估计数应精准等于候选数")

        best_score, best_seq = max(scored)
        self.assertEqual(best_score, max(s for s, _ in scored))
        self.assertGreaterEqual(best_score, -0.7134)
        self.assertTrue(cons.satisfies(best_seq))
        d = obj.evaluate(best_seq)
        self.assertEqual(d["eval_count"], 21)


if __name__ == "__main__":
    unittest.main(verbosity=2)
