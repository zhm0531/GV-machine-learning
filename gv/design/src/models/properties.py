"""真实物理化学性质计算（纯 stdlib，无外部依赖）

全部使用已发表的标准氨基酸量表，不发明权重、不做加权和聚合：
    EISENBERG      Eisenberg et al. 1984, J. Mol. Biol. 共识疏水性量表（归一化）
    CHOU_FASMAN_PB Chou & Fasman 1978, β-折叠倾向
    A3V            AGGRESCAN (Conchillo-Solé et al. 2007) 聚集倾向量表
    KYTE_DOOLITTLE Kyte & Doolittle 1982, 亲水性
    PKA            EMBOSS 默认 pKa，用于等电点

GvpA 是 β-sheet 刚性壳体，内表面疏水、外表面亲水。主目标性质取
signed_amphiphilicity —— 带方向的疏水性矩投影（见该函数说明）。
"""
import math

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

EISENBERG = {
    "A": 0.62, "C": 0.29, "D": -0.90, "E": -0.74, "F": 1.19, "G": 0.48,
    "H": -0.40, "I": 1.38, "K": -1.50, "L": 1.06, "M": 0.64, "N": -0.78,
    "P": 0.12, "Q": -0.85, "R": -2.53, "S": -0.18, "T": -0.05, "V": 1.08,
    "W": 0.81, "Y": 0.26,
}
CHOU_FASMAN_PB = {
    "A": 0.83, "C": 1.19, "D": 0.54, "E": 0.37, "F": 1.38, "G": 0.75,
    "H": 0.87, "I": 1.60, "K": 0.74, "L": 1.30, "M": 1.05, "N": 0.89,
    "P": 0.55, "Q": 1.10, "R": 0.93, "S": 0.75, "T": 1.19, "V": 1.70,
    "W": 1.37, "Y": 1.47,
}
A3V = {
    "A": 0.748, "C": 0.297, "D": -0.860, "E": -0.774, "F": 0.991, "G": 0.097,
    "H": 0.061, "I": 1.161, "K": -0.922, "L": 1.071, "M": 0.728, "N": -0.550,
    "P": -1.115, "Q": -0.598, "R": -0.500, "S": -0.231, "T": -0.143,
    "V": 1.038, "W": 0.708, "Y": 0.714,
}
KYTE_DOOLITTLE = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2,
    "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5,
    "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}
PKA = {"Cterm": 3.6, "Nterm": 8.6, "C": 8.5, "D": 3.9, "E": 4.1,
       "H": 6.5, "K": 10.8, "R": 12.5, "Y": 10.1}
_POS = {"K": 1, "R": 1, "H": 1, "Nterm": 1}
_NEG = {"D": -1, "E": -1, "Cterm": -1, "Y": -1, "C": -1}


def hydrophobic_moment(seq, angle_deg=180.0, window=11):
    """疏水性矩 μH（Eisenberg 1984）。

    μH = (1/n) * |Σ H_i * (cos(iδ), sin(iδ))|
    δ=180° 对应 β-折叠（每残基转 180°），δ=100° 对应 α-螺旋。
    沿序列取滑动窗口后取均值，使局部两亲性变化能反映到分数上。

    返回：各窗口 μH 的均值（≥0，越大越两亲）
    """
    n = len(seq)
    if n == 0:
        return 0.0
    w = min(window, n)
    d = math.radians(angle_deg)
    tot = 0.0
    cnt = 0
    for s in range(n - w + 1):
        x = y = 0.0
        for k in range(w):
            h = EISENBERG[seq[s + k]]
            x += h * math.cos(k * d)
            y += h * math.sin(k * d)
        tot += math.hypot(x, y) / w
        cnt += 1
    return tot / cnt if cnt else 0.0


def beta_propensity(seq):
    """Chou-Fasman β-折叠倾向均值。>1 倾向 β，<1 倾向其他。"""
    return sum(CHOU_FASMAN_PB[a] for a in seq) / len(seq)


def moment_vector(seq, angle_deg=180.0, window=11):
    """各滑动窗口疏水性矩向量的均值（未取向）。用于定义两亲性「方向」。"""
    n = len(seq)
    w = min(window, n)
    d = math.radians(angle_deg)
    sx = sy = 0.0
    cnt = 0
    for s in range(n - w + 1):
        x = y = 0.0
        for k in range(w):
            h = EISENBERG[seq[s + k]]
            x += h * math.cos(k * d)
            y += h * math.sin(k * d)
        sx += x / w
        sy += y / w
        cnt += 1
    return (sx / cnt, sy / cnt) if cnt else (0.0, 0.0)


def signed_amphiphilicity(seq, ref_dir, angle_deg=180.0, window=11):
    """带方向的疏水性矩：把候选的矩向量投影到参考序列的方向上。

    为什么不用幅度 μH：μH = |ΣH·(cos,sin)| 与方向无关，把「内面疏水」换成
    「内面亲水」同样能抬高它。实测参考序列可改区单点最优全部落在 Arg
    (Eisenberg -2.53)，正是把强亲水残基放到了同一相位 —— 幅度上升但两亲性
    方向可能翻转。任务书要求「内表面疏水、外表面亲水」是有方向的，故取投影。

    投影为负 = 两亲性方向与野生型相反，应被惩罚。
    """
    x, y = moment_vector(seq, angle_deg, window)
    rx, ry = ref_dir
    rn = math.hypot(rx, ry)
    if rn == 0:
        return 0.0
    return (x * rx + y * ry) / rn


def aggregation(seq, window=7):
    """AGGRESCAN a3v 聚集倾向。

    返回 (全序列均值, 最强窗口均值)。a3v > 0 偏聚集。
    """
    n = len(seq)
    mean = sum(A3V[a] for a in seq) / n
    w = min(window, n)
    best = max(sum(A3V[seq[s + k]] for k in range(w)) / w
               for s in range(n - w + 1))
    return mean, best


def gravy(seq):
    """Kyte-Doolittle 平均亲水性。>0 疏水，<0 亲水。"""
    return sum(KYTE_DOOLITTLE[a] for a in seq) / len(seq)


def net_charge(seq, ph=7.0):
    """给定 pH 下的净电荷。"""
    q = 0.0
    for a in seq:
        if a in _POS:
            q += 1.0 / (1.0 + 10 ** (ph - PKA[a]))
        elif a in _NEG:
            q -= 1.0 / (1.0 + 10 ** (PKA[a] - ph))
    q += 1.0 / (1.0 + 10 ** (ph - PKA["Nterm"]))
    q -= 1.0 / (1.0 + 10 ** (PKA["Cterm"] - ph))
    return q


def isoelectric_point(seq):
    """等电点 pI：净电荷为 0 的 pH（二分法）。"""
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def all_properties(seq):
    """一次性算全部指标。主目标 target=μH(β)。"""
    agg_mean, agg_max = aggregation(seq)
    return {
        "muH_beta": round(hydrophobic_moment(seq, 180.0), 4),
        "muH_helix": round(hydrophobic_moment(seq, 100.0), 4),
        "beta_propensity": round(beta_propensity(seq), 4),
        "aggregation_mean": round(agg_mean, 4),
        "aggregation_max": round(agg_max, 4),
        "gravy": round(gravy(seq), 4),
        "pI": round(isoelectric_point(seq), 2),
        "net_charge_pH7": round(net_charge(seq, 7.0), 3),
    }


if __name__ == "__main__":
    ref = ("MAVEKVNSSSSLAEVVDRILDKGIVIDAWVRVSLVGIELLAIEARIVIAS"
           "VETYLKYAEAVGLTSQAAVPSA")
    p = all_properties(ref)
    for k, v in p.items():
        print(f"  {k:20} {v}")
    assert p["muH_beta"] > 0, "μH 必须为正"

    # 带符号投影：把序列整体颠倒相位（逆序不改组成）应当掉到负值
    rd = moment_vector(ref)
    assert signed_amphiphilicity(ref, rd) > 0, "参考序列自投影必须为正"
    flipped = ref[1:] + ref[0]
    print(f"\n  参考自投影 {signed_amphiphilicity(ref, rd):.4f}"
          f"  平移一位后 {signed_amphiphilicity(flipped, rd):.4f}")
    print("\n自检通过")
