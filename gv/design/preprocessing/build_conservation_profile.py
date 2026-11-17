"""M2 前置：从 GvpA 天然序列集构建「参考序列锚定的」保守性图谱。

为什么需要这个脚本：
- 包内 `data/GvpA_msa.a3m` 不是比对结果。实测 856 条记录全部长度 60、gap 数 0、
  小写 0，即 create_simple_msa.py 把每条序列截断到最短长度(60)拼出来的假 MSA。
  用它按列统计保守性 = 统计前缀截断伪影，无生物学意义。
- 本机无 mafft/clustalo/muscle/hmmer/biopython（已实测），因此用纯 stdlib 实现
  全局比对(Needleman-Wunsch)，把每条天然序列比对到参考序列，堆叠成参考位点图谱。

产物：data/conservation_profile.json（首次生成约数秒，之后 objective 直接读取）

用法：
    python3 preprocessing/build_conservation_profile.py
    GVP_REF_ID=WP_007357596.1 GVP_MASK_THR=0.9 python3 preprocessing/build_conservation_profile.py
"""
import json
import os
import sys

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
_AA_IDX = {a: i for i, a in enumerate(ALPHABET)}

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("GVP_DATA", os.path.join(HERE, "..", "data"))
OUT = os.path.join(DATA, "conservation_profile.json")

REF_ID = os.environ.get("GVP_REF_ID", "WP_007357596.1")
MASK_THR = float(os.environ.get("GVP_MASK_THR", "0.9"))

MATCH, MISMATCH, GAP = 2, -1, -2


def read_fasta(path):
    """返回 [(id, seq)]，序列大写、仅保留 20 种标准氨基酸。"""
    out, cur, buf = [], None, []
    with open(path) as f:
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
    """全局比对 query 到 ref，返回 (q_aln, r_aln)，长度相等，'-' 表示空位。"""
    n, m = len(query), len(ref)
    if n == 0 or m == 0:
        return query.ljust(m, "-"), ref.ljust(n, "-")
    # traceback: 0=diag, 1=up(gap in ref), 2=left(gap in query)
    tb = [bytearray(m + 1) for _ in range(n + 1)]
    prev = [GAP * j for j in range(m + 1)]
    for j in range(1, m + 1):
        tb[0][j] = 2
    for i in range(1, n + 1):
        cur = [GAP * i] + [0] * m
        tb[i][0] = 1
        qi = query[i - 1]
        for j in range(1, m + 1):
            d = prev[j - 1] + (MATCH if qi == ref[j - 1] else MISMATCH)
            u = prev[j] + GAP
            l = cur[j - 1] + GAP
            if d >= u and d >= l:
                cur[j], tb[i][j] = d, 0
            elif u >= l:
                cur[j], tb[i][j] = u, 1
            else:
                cur[j], tb[i][j] = l, 2
        prev = cur
    # traceback
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


def main():
    fasta = os.path.join(DATA, "GvpA.fasta")
    if not os.path.exists(fasta):
        sys.exit(f"[FATAL] 找不到 {fasta}")
    records = read_fasta(fasta)
    if not records:
        sys.exit(f"[FATAL] {fasta} 无有效序列")

    ref = dict(records).get(REF_ID)
    if ref is None:
        sys.exit(f"[FATAL] 参考序列 {REF_ID} 不在 {fasta} 中；"
                 f"候选 ID: {[h for h, _ in records[:5]]} ...")
    L = len(ref)

    # 每条天然序列比对到参考，堆叠成参考位点列频
    counts = [[0] * len(ALPHABET) for _ in range(L)]
    n_col = [0] * L
    aligned_rows = []
    for _hid, seq in records:
        qaln, raln = nw_align(seq, ref)
        # raln 中非空位的位置就是对齐到参考的位点，按出现顺序回填
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

    mask = [i for i in range(L) if maj_freq[i] >= MASK_THR]

    # 天然性质包络：用 856 条天然序列实测各项性质的取值范围，
    # 作为优化可行域的边界（数据驱动，不设主观阈值）。
    sys.path.insert(0, os.path.join(HERE, "..", "src"))
    from models.properties import (aggregation, beta_propensity, gravy,
                                   isoelectric_point, net_charge)

    envelope = {}
    props = {"aggregation_mean": [], "aggregation_max": [], "gravy": [],
             "beta_propensity": [], "pI": [], "net_charge_pH7": []}
    for row in aligned_rows:
        s = row.replace("-", "")
        if len(s) != L:
            continue
        am, ax = aggregation(s)
        props["aggregation_mean"].append(am)
        props["aggregation_max"].append(ax)
        props["gravy"].append(gravy(s))
        props["beta_propensity"].append(beta_propensity(s))
        props["pI"].append(isoelectric_point(s))
        props["net_charge_pH7"].append(net_charge(s, 7.0))
    for k, v in props.items():
        if not v:
            continue
        v.sort()
        envelope[k] = {"min": round(v[0], 4), "max": round(v[-1], 4),
                       "p5": round(v[int(0.05 * (len(v) - 1))], 4),
                       "p95": round(v[int(0.95 * (len(v) - 1))], 4)}

    payload = {
        "ref_id": REF_ID,
        "ref_seq": ref,
        "n_seq": len(records),
        "alphabet": ALPHABET,
        "profile": profile,
        "majority": "".join(majority),
        "majority_freq": maj_freq,
        "coverage": coverage,
        "mask_threshold": MASK_THR,
        "mask": mask,
        "aligned": aligned_rows,
        "envelope": envelope,
        "wildtype_props": {
            "aggregation_mean": round(aggregation(ref)[0], 4),
            "aggregation_max": round(aggregation(ref)[1], 4),
            "gravy": round(gravy(ref), 4),
            "beta_propensity": round(beta_propensity(ref), 4),
            "pI": round(isoelectric_point(ref), 4),
            "net_charge_pH7": round(net_charge(ref, 7.0), 4),
        },
        "source": "GvpA.fasta, Needleman-Wunsch 参考锚定比对（纯 stdlib）",
    }
    with open(OUT, "w") as f:
        json.dump(payload, f)
    print(f"[OK] {OUT}")
    print(f"  参考 {REF_ID} len={L}  序列数={len(records)}")
    print(f"  平均列覆盖={sum(coverage)/L:.3f}  平均多数残基频率={sum(maj_freq)/L:.3f}")
    print(f"  锁定(硬约束候选, freq>={MASK_THR}) {len(mask)}/{L} 位: {mask}")
    print(f"  可改位点 {L - len(mask)} 个")
    print(f"  天然性质包络:")
    for k, v in envelope.items():
        print(f"    {k:18} [{v['min']:>8.4f}, {v['max']:>8.4f}]"
              f"   P5={v['p5']:>8.4f} P95={v['p95']:>8.4f}")


if __name__ == "__main__":
    main()
