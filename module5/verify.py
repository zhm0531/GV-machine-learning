"""模块5 · 块1：约束验证（Constraint Verification）。

验证优化器产出的候选序列是否仍满足 GvpA 硬约束。只读序列、算指标、出结果。

两级判定：
  1. fast_check()：纯 Python，微秒级
     - 长度 == 参考长度（substitution-only）
     - 仅含 20 种标准氨基酸
     - 掩码(MASK)锁定位点与参考完全一致
  2. heavy_check()：HMMER 调 PF00741 --cut_ga，秒级，只跑通过快速判定的 top-k

指标：
  - 约束满足率 R  = 通过判定的序列数 / 总序列数
  - 保守位点保持率 C = 1 - 锁定改动总次数 / (序列数 × 锁定数)

输出：verification.csv + constraint_report.txt
"""
import argparse
import os
import sys

from . import common

DEFAULT_INPUT = "best.fasta"
DEFAULT_OUTDIR = "verify_out"
DEFAULT_MASK = "mask.txt"
DEFAULT_REF = "ref.txt"
DEFAULT_TOP_K = 50


def _load_pfam_scan():
    """尝试从 gv 仓库加载 pfam_scan，失败返回 None。"""
    mod = common.load_optional_module("models.constraints", ["pfam_scan"])
    return getattr(mod, "pfam_scan", None) if mod else None


def _load_ref_seq():
    """尝试从 gv 仓库加载参考序列，失败返回 None。"""
    mod = common.load_optional_module("models.constraints", ["REF_SEQ"])
    return getattr(mod, "REF_SEQ", None) if mod else None


def fast_check(seq, ref_seq, mask):
    """快速判定，返回违规描述列表；空列表 = 通过。"""
    violations = []
    if not isinstance(seq, str):
        return ["非字符串输入: %s" % type(seq).__name__]
    if len(seq) != len(ref_seq):
        violations.append("长度 %d != 参考长度 %d" % (len(seq), len(ref_seq)))
    bad = set(seq) - set(common.ALPHABET)
    if bad:
        violations.append("非法残基 %s" % sorted(bad))
    if len(seq) == len(ref_seq) and not bad and mask:
        diff = [i for i in mask if seq[i] != ref_seq[i]]
        if diff:
            violations.append("锁定位点被改动 %d 处: %s" % (len(diff), diff[:10]))
    return violations


def heavy_check(seqs, names):
    """重判定：HMMER pfam_scan。未接入或失败时逐条返回 None。"""
    pfam_scan = _load_pfam_scan()
    if pfam_scan is None:
        print("[警告] 未接入 gv 的 pfam_scan，重判定跳过", file=sys.stderr)
        return [None] * len(seqs)
    try:
        return pfam_scan(seqs, names=names)
    except Exception as exc:
        print("[警告] pfam_scan 执行失败: %s" % exc, file=sys.stderr)
        return [None] * len(seqs)


def compute_metrics(records, ref_seq, mask, top_k=DEFAULT_TOP_K):
    """计算逐条明细与汇总指标。"""
    rows = []
    mask_set = set(mask)
    for i, (name, seq) in enumerate(records):
        violations = fast_check(seq, ref_seq, mask)
        passed_fast = len(violations) == 0
        mutated_locked = (
            sum(1 for m in mask_set if seq[m] != ref_seq[m])
            if len(seq) == len(ref_seq) else None)
        rows.append({
            "rank": i + 1,
            "name": name,
            "length": len(seq),
            "passed_fast": passed_fast,
            "mutated_locked_sites": mutated_locked,
            "violations": "; ".join(violations),
            "passed_pfam": None,
            "pfam_evalue": None,
            "pfam_score": None,
        })

    cand = [(r["name"], records[i][1]) for i, r in enumerate(rows)
            if r["passed_fast"]][:top_k]
    if cand:
        names = [c[0] for c in cand]
        seqs = [c[1] for c in cand]
        results = heavy_check(seqs, names)
        name2res = dict(zip(names, results))
        for r in rows:
            res = name2res.get(r["name"])
            if res is not None:
                r["passed_pfam"] = res.get("passed")
                r["pfam_evalue"] = res.get("evalue")
                r["pfam_score"] = res.get("score")

    n = len(rows)
    pfam_available = any(r["passed_pfam"] is not None for r in rows)
    if n == 0:
        R = 0.0
    elif pfam_available:
        R = sum(1 for r in rows if r["passed_fast"] and r["passed_pfam"]) / n
    else:
        R = sum(1 for r in rows if r["passed_fast"]) / n

    locked_total = len(mask_set)
    total_mut = sum((r["mutated_locked_sites"] or 0) for r in rows)
    total_locked_checks = locked_total * n
    C = 1.0 - (total_mut / total_locked_checks) if total_locked_checks else 1.0

    metrics = {
        "n_sequences": n,
        "constraint_satisfaction_R": round(R, 4),
        "locked_sites": locked_total,
        "seqs_with_locked_mutation": sum(1 for r in rows if r["mutated_locked_sites"]),
        "locked_site_mutations_total": total_mut,
        "retention_C": round(C, 4),
        "fast_pass_count": sum(1 for r in rows if r["passed_fast"]),
        "pfam_pass_count": sum(1 for r in rows if r["passed_pfam"]),
        "pfam_available": pfam_available,
    }
    return rows, metrics


def write_outputs(rows, metrics, outdir=DEFAULT_OUTDIR):
    """写出 verification.csv 与 constraint_report.txt。"""
    common.write_csv(os.path.join(outdir, "verification.csv"), rows)
    m = metrics
    lines = [
        "=" * 60,
        "约束验证报告 (模块5 - 块1)",
        "=" * 60,
        "序列总数            : %d" % m["n_sequences"],
        "快速判定通过         : %d" % m["fast_pass_count"],
        "Pfam(PF00741) 通过   : %d%s" % (
            m["pfam_pass_count"], "" if m["pfam_available"] else "  (重判定未执行)"),
        "约束满足率 R         : %.2f%%" % (m["constraint_satisfaction_R"] * 100),
        "锁定位点数           : %d" % m["locked_sites"],
        "含锁定改动的序列数   : %d" % m["seqs_with_locked_mutation"],
        "锁定改动总次数       : %d" % m["locked_site_mutations_total"],
        "保守位点保持率 C     : %.2f%%" % (m["retention_C"] * 100),
        "-" * 60,
        "逐条明细见 verification.csv",
        "=" * 60,
    ]
    common.write_report(os.path.join(outdir, "constraint_report.txt"), lines)
    print("\n".join(lines))
    print("\n明细已写入: %s" % os.path.join(outdir, "verification.csv"))


def resolve_ref(args):
    """确定参考序列：优先 --ref 文件，其次 gv 仓库，最后给出占位并警告。"""
    if os.path.exists(args.ref):
        return "".join(common.iter_sequence_lines(args.ref)).split()[0]
    ref_seq = _load_ref_seq()
    if ref_seq is None:
        print("[提示] 未提供参考序列，使用占位（结果不可用于结论）", file=sys.stderr)
        return "M" * 72
    return ref_seq


def main(argv=None):
    ap = argparse.ArgumentParser(description="约束验证（模块5 - 块1）")
    ap.add_argument("--input", default=DEFAULT_INPUT, help="候选序列 fasta")
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR, help="输出目录")
    ap.add_argument("--mask", default=DEFAULT_MASK, help="掩码文件（0 基索引）")
    ap.add_argument("--ref", default=DEFAULT_REF, help="参考序列文件")
    ap.add_argument("--top-k", type=int, default=DEFAULT_TOP_K,
                    help="重判定只跑前 k 条")
    args = ap.parse_args(argv)

    common.bootstrap_src()
    ref_seq = resolve_ref(args)
    # 掩码：显式 --mask 优先，否则用项目 mod1 正式掩码
    mask_path = (args.mask if args.mask and os.path.exists(args.mask)
                 else common.ensure_mask_path())
    mask = common.read_locked_positions(mask_path) if mask_path else []
    records = common.read_fasta(args.input)
    if not records:
        print("未读到任何序列，退出", file=sys.stderr)
        return
    rows, metrics = compute_metrics(records, ref_seq, mask, top_k=args.top_k)
    write_outputs(rows, metrics, outdir=args.outdir)


if __name__ == "__main__":
    main()