# -*- coding: utf-8 -*-
"""一键运行模块5 全流程：约束验证 -> 副作用分析 -> 中期报告。

自动定位项目根目录下的正式掩码（module1）、数据目录（gv/design/data）和
接口文件（module2），无需再手动指定路径；也保留显式覆盖参数。

用法：
    python run_all.py --input best.fasta
    python run_all.py --input best.fasta --pfam
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from module5 import verify, side_effect, report
from module5 import common


def main():
    root = common.PROJECT_ROOT
    default_data = common.find_project_data() or os.path.join(root, "gv", "design", "data")
    default_mask = common.ensure_mask_path() or "mask.txt"

    ap = argparse.ArgumentParser(description="模块5 全流程（自动定位路径）")
    ap.add_argument("--input", default="best.fasta")
    ap.add_argument("--data", default=default_data)
    ap.add_argument("--mask", default=default_mask)
    ap.add_argument("--ref", default="")
    ap.add_argument("--history", default="history.csv",
                    help="模块4 history.csv 或目录（可选）")
    ap.add_argument("--pfam", action="store_true",
                    help="对通过快速判定的候选执行 HMMER PF00741 重判")
    ap.add_argument("--verify-out", default="verify_out")
    ap.add_argument("--side-out", default="sideeffect_out")
    af = ap.parse_args()

    if not af.data:
        print("错误：未找到项目数据目录（conservation_profile.json），请用 --data 指定。",
              file=sys.stderr)
        return 1

    print("[1/3] 约束验证")
    verify.main(["--input", af.input, "--mask", af.mask,
                 "--outdir", af.verify_out] +
                (["--ref", af.ref] if af.ref else []))
    print("\n[2/3] 副作用分析")
    side_effect.main(["--input", af.input, "--data", af.data,
                      "--mask", af.mask, "--outdir", af.side_out])
    print("\n[3/3] 中期报告")
    report.main(["--verify", os.path.join(af.verify_out, "verification.csv"),
                 "--side", os.path.join(af.side_out, "side_effect.csv"),
                 "--history", af.history])
    print("\n全部完成。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
