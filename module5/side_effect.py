"""模块5 · 块2：副作用分析（E4），回答 RQ3：优化是否引入非预期变化。

逐条候选检查：
  1. 位点审计：改动落在可变区还是非可变区；非可变区改动比例 >10% 报警
  2. 锁定检查：掩码位点必须 0 改动
  3. 接口检查：GvpC 组装接口位点改动必须 0（探索项）
  4. 理化性质偏移：与参考对比，是否越出天然包络 C_env
  5. 结构稳定性代理：折叠代理分相对参考的下降（>0.05 报警）
     （无 AlphaFold/ColabFold 时用代理；接入真实 pLDDT 时替换 fold_proxy 即可）

输出：side_effect.csv + side_effect_report.txt
"""
import argparse
import json
import os

from . import common
import csv

DEFAULT_INPUT = "best.fasta"
DEFAULT_DATA = ""
DEFAULT_PROFILE = "conservation_profile.json"
DEFAULT_MASK = "mask.txt"
DEFAULT_IFACE = "domain_interface_sites.csv"
DEFAULT_OUTDIR = "sideeffect_out"

VARIABLE_THR = 0.70
OUTSIDE_RATIO_LIMIT = 0.10
FOLD_DROP_LIMIT = 0.05
PROP_KEYS = ["aggregation_mean", "aggregation_max", "gravy",
             "beta_propensity", "pI", "net_charge_pH7"]


def _properties_module():
    """加载 properties 模块（含 all_properties），失败返回 None。"""
    return common.load_optional_module("models.properties", ["all_properties"])


def _safe_all_properties(seq):
    """调用 properties.all_properties；模块缺失时抛 RuntimeError（带清晰提示）。"""
    mod = _properties_module()
    if mod is None or not hasattr(mod, "all_properties"):
        raise RuntimeError(
            "找不到 properties 模块（无法计算理化性质）。请设置 GVP_SRC "
            "指向 gv/design/src，或在 --data 提供 GvpA 数据。")
    return mod.all_properties(seq)


def _resolve_data_dir(data, here):
    """确定数据目录：优先显式 --data；否则项目根自动发现；最后工作目录 data/。"""
    src_dir = common.bootstrap_src()
    candidates = [data,
                  common.find_project_data(),
                  (os.path.join(src_dir, '..', 'data') if src_dir else ''),
                  os.path.join(here, 'data')]
    for cand in candidates:
        if cand and os.path.isfile(os.path.join(cand, DEFAULT_PROFILE)):
            return cand
    raise FileNotFoundError(
        '找不到数据目录（需含 %s）。请用 --data 指定。' % DEFAULT_PROFILE)


def load_interface_sites(path):
    """读 GvpC 接口位点：支持三种来源。

    1. JSON 对象 {"sites": [...]}（兼容旧格式）
    2. 纯文本索引行（0 基，支持 # 注释与 start-end 范围）
    3. module2 的 domain_interface_sites.csv（列 start/end，按 1 基位置取整段位点）
    不存在返回 None。
    """
    if not path or not os.path.exists(path):
        return None
    if path.lower().endswith('.csv'):
        return _load_iface_csv(path)
    with open(path, 'r', encoding='utf-8') as f:
        raw = f.read().strip()
    if not raw:
        return None
    if raw.startswith('{'):
        return [int(x) for x in json.loads(raw).get('sites', [])]
    return _parse_iface_index_lines(raw)


def _parse_iface_index_lines(raw):
    sites = set()
    for ln in raw.splitlines():
        ln = ln.split('#')[0].strip()
        if not ln:
            continue
        if '-' in ln:
            try:
                a, b = ln.split('-')
                sites.update(range(int(a), int(b) + 1))
            except ValueError:
                pass
        else:
            try:
                sites.add(int(ln))
            except ValueError:
                pass
    return sorted(sites)


def _load_iface_csv(path):
    """读 module2 接口 CSV：取 type 含 assembly/interface/domain 的行，start/end 转 0 基位点。"""
    sites = set()
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            st = row.get('start'); en = row.get('end')
            typ = (row.get('type') or '').lower()
            if st is None or en is None or st == '' or en == '':
                continue
            try:
                a, b = int(st), int(en)
            except ValueError:
                continue
            # PF00741 domain / assembly interface 都纳入；prior_work 仅参考不纳入
            if 'prior' in typ:
                continue
            sites.update(range(a - 1, b))  # CSV 是 1 基，转 0 基
    return sorted(sites) if sites else None




def fold_proxy(seq, ref_props):
    """结构稳定性代理：关键结构性质相对参考的接近度，1=与参考一致。"""
    p = _safe_all_properties(seq)
    keys = ["beta_propensity", "aggregation_mean", "gravy", "pI"]
    dev = sum(abs(p[k] - ref_props[k]) / (abs(ref_props[k]) + 1e-6) for k in keys)
    return max(0.0, 1.0 - dev / len(keys))


def analyze(records, prof, mask, iface):
    """对候选序列逐条做副作用分析，返回明细行列表。"""
    ref = prof["ref_seq"]
    env = prof.get("envelope", {})
    maj_freq = prof.get("majority_freq", [1.0] * len(ref))
    variable = set(i for i in range(len(ref)) if maj_freq[i] <= VARIABLE_THR)
    nonvar = set(range(len(ref))) - variable
    mask_set = set(mask)
    ref_props = _safe_all_properties(ref)
    ref_fold = 1.0

    rows = []
    for name, seq in records:
        r = {"name": name, "length": len(seq)}
        if len(seq) != len(ref):
            r.update({"verdict": "长度异常", "n_edits": None})
            rows.append(r)
            continue

        edits = [i for i in range(len(ref)) if seq[i] != ref[i]]
        edits_in_var = [i for i in edits if i in variable]
        edits_out_var = [i for i in edits if i in nonvar]
        outside_ratio = len(edits_out_var) / len(edits) if edits else 0.0
        locked_changes = [i for i in edits if i in mask_set]
        iface_changes = [i for i in edits if i in set(iface)] if iface else []

        p = _safe_all_properties(seq)
        prop_flags = []
        for k in PROP_KEYS:
            e = env.get(k)
            if e and not (e["min"] <= p[k] <= e["max"]):
                prop_flags.append("%s越界" % k)
        fp = fold_proxy(seq, ref_props)
        fold_drop = ref_fold - fp

        flags = []
        if locked_changes:
            flags.append("锁定改动%d" % len(locked_changes))
        if outside_ratio > OUTSIDE_RATIO_LIMIT:
            flags.append("非可变区改动%.0f%%" % (outside_ratio * 100))
        if iface and iface_changes:
            flags.append("接口改动%d" % len(iface_changes))
        if fold_drop > FOLD_DROP_LIMIT:
            flags.append("折叠代理降%.3f" % fold_drop)
        if prop_flags:
            flags.append(",".join(prop_flags))

        r.update({
            "n_edits": len(edits),
            "edits_in_variable": len(edits_in_var),
            "edits_outside_variable": len(edits_out_var),
            "outside_ratio": round(outside_ratio, 3),
            "locked_changes": len(locked_changes),
            "interface_changes": (len(iface_changes) if iface else None),
            "fold_proxy": round(fp, 4),
            "fold_drop": round(fold_drop, 4),
            "prop_flags": ",".join(prop_flags),
            "verdict": ("OK" if not flags else " | ".join(flags)),
        })
        rows.append(r)
    return rows


def write_outputs(rows, outdir, iface):
    """写出 side_effect.csv 与 side_effect_report.txt。"""
    common.write_csv(os.path.join(outdir, "side_effect.csv"), rows)
    n = len(rows)
    abnormal = [r for r in rows if str(r.get("verdict")) not in ("OK", "长度异常")]
    lines = [
        "=" * 60,
        "副作用分析报告 (模块5 - 块2)",
        "=" * 60,
        "候选序列数        : %d" % n,
        "存在副作用的序列   : %d" % len(abnormal),
        "非可变区改动阈值   : >%.0f%%" % (OUTSIDE_RATIO_LIMIT * 100),
        "折叠代理下降阈值   : >%.3f" % FOLD_DROP_LIMIT,
        "GvpC 接口位点      : %s" % ("未提供" if iface is None else "%d 个" % len(iface)),
        "-" * 60,
    ]
    for r in rows:
        lines.append("  %-20s edits=%s 非可变区=%s 锁定改动=%s 接口=%s 折叠降=%s  %s" % (
            r.get("name"), r.get("n_edits"), r.get("edits_outside_variable"),
            r.get("locked_changes"), r.get("interface_changes"),
            r.get("fold_drop"), r.get("verdict")))
    lines += ["=" * 60, "明细见 side_effect.csv", "=" * 60]
    common.write_report(os.path.join(outdir, "side_effect_report.txt"), lines)
    print("\n".join(lines))





def _default_mask_path():
    """返回本项目 module1 正式掩码路径，存在则返回，否则返回 ''。"""
    for cand in (os.path.join(common.PROJECT_ROOT, "module1", "outputs", "mask.txt"),
                 os.path.join(common.PROJECT_ROOT, "gv", "design", "src", "module4", "config", "official_mask.txt")):
        if os.path.isfile(cand):
            return cand
    return ""


def resolve_iface_path(iface, data_dir):
    """解析接口文件路径：显式传参优先；否则在 module2 与数据目录内查找。"""
    if iface and os.path.exists(iface):
        return iface
    for cand in (os.path.join(common.PROJECT_ROOT, "module2", "out", "domain_interface_sites.csv"),
                 os.path.join(common.PROJECT_ROOT, "module2", "out", "gvpc_interface_sites.json"),
                 os.path.join(data_dir, iface) if iface else "",
                 os.path.join(data_dir, "domain_interface_sites.csv")):
        if cand and os.path.exists(cand):
            return cand
    return iface


def main(argv=None):
    ap = argparse.ArgumentParser(description="副作用分析（模块5 - 块2）")
    ap.add_argument("--input", default=DEFAULT_INPUT)
    ap.add_argument("--data", default=DEFAULT_DATA)
    ap.add_argument("--profile", default=DEFAULT_PROFILE)
    ap.add_argument("--mask", default=DEFAULT_MASK)
    ap.add_argument("--iface", default=DEFAULT_IFACE)
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR)
    args = ap.parse_args(argv)

    here = os.path.dirname(os.path.abspath(__file__))
    data_dir = _resolve_data_dir(args.data, here)
    prof = common.load_profile(data_dir, args.profile)
    # 掩码优先级：显式 --mask > 项目 module1 正式掩码 > profile.mask
    mask_path = (args.mask if args.mask and os.path.exists(args.mask)
                 else common.ensure_mask_path())
    mask = (common.read_locked_positions(mask_path)
            if mask_path else prof.get("mask", []))
    # 接口文件若显式指定则用；否则在本项目 module2 与数据目录内查找 CSV/JSON
    iface = load_interface_sites(resolve_iface_path(args.iface, data_dir))
    records = common.read_fasta(args.input)
    rows = analyze(records, prof, mask, iface)
    write_outputs(rows, args.outdir, iface)


if __name__ == "__main__":
    main()