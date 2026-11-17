# -*- coding: utf-8 -*-
"""
用法：
    python scripts/build_profile.py
    python scripts/build_profile.py --ref-id WP_007357596.1 --mask-thr 0.90
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "outputs")

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
_AA_IDX = {a: i for i, a in enumerate(ALPHABET)}

DEFAULT_REF_ID = os.environ.get("GVP_REF_ID", "WP_007357596.1")
DEFAULT_MASK_THR = float(os.environ.get("GVP_MASK_THR", "0.90"))

MATCH, MISMATCH, GAP = 2, -1, -2

def read_fasta(path):
    """返回 [(id, seq)]；序列大写、仅保留 20 种标准氨基酸。"""
    out, cur, buf = [], None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    out.append((cur, "".join(buf)))
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append("".join(c for c in line.upper() if c in _AA_IDX))
    if cur:
        out.append((cur, "".join(buf)))
    return out

def nw_align(query, ref):
    """全局比对 query 到 ref，返回 (q_aln, r_aln)；等长，- 为 indel。"""
    n, m = len(query), len(ref)
    if n == 0 or m == 0:
        return query.ljust(m, "-"), ref.ljust(n, "-")
    tb = [bytearray(m + 1) for _ in range(n + 1)]
    prev = [GAP * j for j in range(m + 1)]
    for j in range(1, m + 1):
        tb[0][j] = 2
    for i in range(1, n + 1):
        cur = [GAP * i] + [0] * m
        tb[i][0] = 1
        qi = query[i - 1]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (MATCH if qi == ref[j - 1] else MISMATCH)
            up = prev[j] + GAP
            left = cur[j - 1] + GAP
            if diag >= up and diag >= left:
                cur[j], tb[i][j] = diag, 0
            elif up >= left:
                cur[j], tb[i][j] = up, 1
            else:
                cur[j], tb[i][j] = left, 2
        prev = cur
    a, b = [], []
    i, j = n, m
    while i > 0 or j > 0:
        t = tb[i][j]
        if i > 0 and j > 0 and t == 0:
            a.append(query[i - 1]); b.append(ref[j - 1]); i -= 1; j -= 1
        elif i > 0 and (j == 0 or t == 1):
            a.append(query[i - 1]); b.append("-"); i -= 1
        else:
            a.append("-"); b.append(ref[j - 1]); j -= 1
    return "".join(reversed(a)), "".join(reversed(b))

def build_profile(records, ref_seq):
    """返回 (profile, majority, maj_freq, coverage, aligned_rows)。"""
    L = len(ref_seq)
    counts = [[0] * len(ALPHABET) for _ in range(L)]
    n_col = [0] * L
    aligned_rows = []
    for _hid, seq in records:
        qaln, raln = nw_align(seq, ref_seq)
        i = 0
        for c, rch in enumerate(raln):
            if rch == "-":
                continue
            q = qaln[c]
            if q in _AA_IDX:
                counts[i][_AA_IDX[q]] += 1
            n_col[i] += 1
            i += 1
        aligned_rows.append(qaln)
    profile, majority, maj_freq, coverage = [], [], [], []
    for i in range(L):
        tot = n_col[i] or 1
        freq = [round(c / tot, 4) for c in counts[i]]
        profile.append(freq)
        k = max(range(len(ALPHABET)), key=lambda x: freq[x])
        majority.append(ALPHABET[k])
        maj_freq.append(freq[k])
        coverage.append(round(n_col[i] / len(records), 3))
    return profile, majority, maj_freq, coverage, aligned_rows


def _hmmer_bin():
    """返回 hmmsearch 可执行路径；找不到返回 None（不抛错）。"""
    exe = os.environ.get("HMMSEARCH_BIN") or shutil.which("hmmsearch")
    return exe or None


def pfam_verify(ref_seq, mask, hmm_path=None):
    """可选加强校验：用 PF00741 + --cut_ga 复核参考序列仍属 GvpA 家族。

    返回 dict：pfam_available / passed / evalue / score / note。
    hmmsearch 或模型缺失时优雅降级（pfam_available=False），不影响主流程。
    """
    hmm_path = hmm_path or os.path.join(DATA, "PF00741.hmm")
    if _hmmer_bin() is None or not os.path.exists(hmm_path):
        return {"pfam_available": False, "passed": None, "evalue": None,
                "score": None,
                "note": "未安装 hmmsearch 或缺少 PF00741.hmm，PF00741 校验跳过"}
    with tempfile.TemporaryDirectory() as td:
        fasta = os.path.join(td, "in.fa")
        out = os.path.join(td, "out.domtbl")
        with open(fasta, "w", encoding="utf-8") as f:
            f.write(">ref\n%s\n" % ref_seq)
        try:
            subprocess.run([_hmmer_bin(), "--noali", "--cut_ga",
                            "--domtblout", out, hmm_path, fasta],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           check=True)
        except Exception as e:
            return {"pfam_available": False, "passed": None, "evalue": None,
                    "score": None, "note": "PF00741 扫描执行失败: %s" % e}
        hits = {}
        with open(out, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                p = line.split()
                if len(p) < 8:
                    continue
                hits.setdefault(p[0], (float(p[6]), float(p[7])))
        hit = hits.get("ref")
        return {"pfam_available": True,
                "passed": hit is not None,
                "evalue": (hit[0] if hit else None),
                "score": (hit[1] if hit else None),
                "note": ""}


def main():
    ap = argparse.ArgumentParser(description="模块1：GvpA 数据准备 + 正式掩码")
    ap.add_argument("--ref-id", default=DEFAULT_REF_ID)
    ap.add_argument("--mask-thr", type=float, default=DEFAULT_MASK_THR)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    fasta = os.path.join(DATA, "GvpA.fasta")
    if not os.path.exists(fasta):
        sys.exit("[FATAL] 找不到 %s" % fasta)

    records = read_fasta(fasta)
    ref_seq = dict(records).get(args.ref_id)
    if ref_seq is None:
        sys.exit("[FATAL] 参考 %s 不在 GvpA.fasta 中" % args.ref_id)
    L = len(ref_seq)

    profile, majority, maj_freq, coverage, aligned_rows = build_profile(records, ref_seq)
    mask = [i for i in range(L) if maj_freq[i] >= args.mask_thr]
    fr_mask = [i for i in range(L) if i not in set(mask)]
    pfam_res = pfam_verify(ref_seq, mask)

    payload = {
        "ref_id": args.ref_id,
        "ref_seq": ref_seq,
        "n_seq": len(records),
        "L": L,
        "alphabet": ALPHABET,
        "profile": profile,
        "majority": "".join(majority),
        "majority_freq": maj_freq,
        "coverage": coverage,
        "mask_threshold": args.mask_thr,
        "mask": mask,
        "free_positions": fr_mask,
        "mask_summary": {
            "locked": len(mask),
            "free": len(fr_mask),
            "mean_majority_freq": round(sum(maj_freq) / L, 4),
            "mean_coverage": round(sum(coverage) / L, 4),
        },
        "aligned": aligned_rows,
        "pfam_verify": pfam_res,
        "source": "GvpA.fasta + 参考锚定 Needleman-Wunsch 比对（纯 stdlib）",
    }
    prof_path = os.path.join(OUT, "conservation_profile.json")
    with open(prof_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    mask_path = os.path.join(OUT, "mask.txt")
    with open(mask_path, "w", encoding="utf-8") as f:
        f.write("# 正式锁定位点掩码 M（0 基索引，多数残基频率 >= %.2f）\n" % args.mask_thr)
        f.write("# ref_id=%s  L=%d  locked=%d  free=%d\n" % (args.ref_id, L, len(mask), len(fr_mask)))
        f.write("\n".join(str(i) for i in mask) + "\n")

    print("[OK] %s" % prof_path)
    print("[OK] %s" % mask_path)
    print("  参考 %s  L=%d  序列数=%d" % (args.ref_id, L, len(records)))
    print("  平均列覆盖=%.3f  平均多数频率=%.3f"
          % (payload["mask_summary"]["mean_coverage"],
             payload["mask_summary"]["mean_majority_freq"]))
    print("  正式掩码(锁定)=%d/%d  可改=%d  阈值>=%.2f"
          % (len(mask), L, len(fr_mask), args.mask_thr))
    print("  锁定位置: %s" % ",".join(map(str, mask)))
    print("  可改位置: %s" % ",".join(map(str, fr_mask)))
    if pfam_res["pfam_available"]:
        print("  校验 PF00741 --cut_ga: passed=%s  E=%.2e  score=%.1f" % (
              pfam_res["passed"], pfam_res["evalue"], pfam_res["score"]))
    else:
        print("  校验 PF00741: %s" % pfam_res["note"])

if __name__ == "__main__":
    main()