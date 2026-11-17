"""模块5 · 块3：中期演示材料生成（M3 收尾）。

汇总模块1-4 的产出，生成中期报告与演示素材：
  - 模块4 history.csv（best-so-far 收敛曲线）-> 汇总各策略/预算的最终值
  - verify 的 verification.csv（约束验证 R/C）
  - side_effect 的 side_effect.csv（副作用）
  - 生成 midterm_report.md；有 matplotlib 时额外输出收敛曲线 PNG

不依赖第三方库即可运行。
"""
import argparse
import csv
import glob
import os
from collections import defaultdict

DEFAULT_HISTORY = "history.csv"
DEFAULT_VERIFY = "verify_out/verification.csv"
DEFAULT_SIDE = "sideeffect_out/side_effect.csv"
DEFAULT_OUT = "midterm_report.md"


def read_csv(path):
    """读 CSV；文件不存在返回空列表。"""
    if not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def find_cols(rows, *cands):
    """在行字段中按候选名（忽略大小写/子串）查找列名。"""
    if not rows:
        return None
    cols = list(rows[0].keys())
    for cand in cands:
        for col in cols:
            if col.lower() == cand.lower():
                return col
    for cand in cands:
        for col in cols:
            if cand.lower() in col.lower():
                return col
    return None


def summarize_history(rows):
    """把 history.csv 汇总成 {(strategy, budget, seed): 最终 best}。"""
    if not rows:
        return {}, None
    c_best = find_cols(rows, "best", "best_so_far", "best_f", "score")
    c_strat = find_cols(rows, "strategy", "optimizer", "algo", "model")
    c_budget = find_cols(rows, "budget", "max_evals")
    c_seed = find_cols(rows, "seed")
    grouped = defaultdict(list)
    for r in rows:
        key = tuple(r.get(k) for k in (c_strat, c_budget, c_seed) if k)
        try:
            grouped[key].append(float(r[c_best]))
        except (TypeError, ValueError, KeyError):
            pass
    return {k: max(v) for k, v in grouped.items() if v}, c_best


def summarize_verify(rows):
    """统计约束验证。"""
    if not rows:
        return {}
    return {
        "n": len(rows),
        "fast_pass": sum(1 for r in rows if str(r.get("passed_fast")).lower() == "true"),
        "pfam_pass": sum(1 for r in rows if str(r.get("passed_pfam")).lower() == "true"),
    }


def summarize_side(rows):
    """统计副作用。"""
    if not rows:
        return {}
    n = len(rows)
    ok = sum(1 for r in rows if str(r.get("verdict")) == "OK")
    return {"n": n, "ok": ok, "abnormal": n - ok}


def maybe_plot(rows, outdir):
    """有 matplotlib 时画收敛曲线；否则返回 None。"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    c_eval = find_cols(rows, "eval", "step", "n", "iteration")
    c_best = find_cols(rows, "best", "best_so_far", "best_f", "score")
    c_strat = find_cols(rows, "strategy", "optimizer", "algo", "model")
    if not (c_eval and c_best):
        return None
    series = defaultdict(list)
    for r in rows:
        try:
            series[r.get(c_strat, "all")].append((float(r[c_eval]), float(r[c_best])))
        except (TypeError, ValueError, KeyError):
            pass
    if not series:
        return None
    plt.figure(figsize=(8, 5))
    for name, pts in series.items():
        pts.sort()
        plt.plot([p[0] for p in pts], [p[1] for p in pts], label=str(name))
    plt.xlabel("evaluations"); plt.ylabel("best-so-far")
    plt.title("Optimization convergence"); plt.legend()
    png = os.path.join(outdir, "convergence.png")
    plt.tight_layout(); plt.savefig(png, dpi=150); plt.close()
    return png


def render(hist_rows, verify_rows, side_rows, out):
    """汇总并写出 midterm_report.md，返回 (报告路径, 曲线路径|None)。"""
    hist_summary, _ = summarize_history(hist_rows)
    v = summarize_verify(verify_rows)
    s = summarize_side(side_rows)
    outdir = os.path.dirname(os.path.abspath(out)) or "."
    png = maybe_plot(hist_rows, outdir) if hist_rows else None

    lines = ["# GV02-04 中期报告（M3）", "", "## 1. 优化收敛（模块4）", ""]
    if hist_summary:
        lines.append("| 策略/预算/seed | 最终 best-so-far |")
        lines.append("|---|---|")
        for key, val in sorted(hist_summary.items(), key=lambda x: -x[1]):
            lines.append("| %s | %.4f |" % (" / ".join(str(x) for x in key), val))
    else:
        lines.append("_未找到 history.csv（模块4 产出）_")
    lines += ["", "## 2. 约束验证（代码1）", "",
              "- 序列数: %s" % v.get("n", "-"),
              "- 快速判定通过: %s" % v.get("fast_pass", "-"),
              "- Pfam(PF00741) 通过: %s" % v.get("pfam_pass", "-"),
              "", "## 3. 副作用分析（代码2）", "",
              "- 序列数: %s" % s.get("n", "-"),
              "- 无副作用(OK): %s" % s.get("ok", "-"),
              "- 存在副作用: %s" % s.get("abnormal", "-"), ""]
    if png:
        lines += ["## 4. 收敛曲线", "", "![convergence](%s)" % os.path.basename(png), ""]
    lines += ["---", "_本报告由 模块5 自动生成；数据来源为模块1-4 的产出文件。_"]
    os.makedirs(outdir, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return out, png


def main(argv=None):
    ap = argparse.ArgumentParser(description="中期演示材料生成（模块5 - 块3）")
    ap.add_argument("--history", default=DEFAULT_HISTORY,
                    help="模块4 history.csv 或目录")
    ap.add_argument("--verify", default=DEFAULT_VERIFY)
    ap.add_argument("--side", default=DEFAULT_SIDE)
    ap.add_argument("--out", default=DEFAULT_OUT)
    args = ap.parse_args(argv)

    hist_rows = []
    if os.path.isdir(args.history):
        for p in glob.glob(os.path.join(args.history, "**", "history.csv"),
                           recursive=True):
            hist_rows += read_csv(p)
    else:
        hist_rows = read_csv(args.history)
    verify_rows = read_csv(args.verify)
    side_rows = read_csv(args.side)

    out, png = render(hist_rows, verify_rows, side_rows, args.out)
    print("已生成: %s" % out)
    if png:
        print("已生成: %s" % png)


if __name__ == "__main__":
    main()