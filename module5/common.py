"""模块5 · 共享工具（common）。

被 verify / side_effect / report 三个子模块复用，消除重复实现：
  - fasta / 掩码 / 序列文件读写
  - 外部 gv 设计仓库（GVP_SRC）的定位与 sys.path 注入
  - 项目数据目录（GVP_DATA）与参考/掩码/包络的加载
  - 统一的文件存在性错误提示

纯标准库实现，无第三方依赖。
"""
import csv
import json
import os
import sys

# --------------------------------------------------------------------------- #
# 常量与默认路径
# --------------------------------------------------------------------------- #
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def find_project_src() -> str:
    """在本项目根目录下定位含 models/ 的 gv/design/src（找不到返回 ''）。"""
    for cand in (
        os.path.join(PROJECT_ROOT, "gv", "design", "src"),
        os.path.join(PROJECT_ROOT, "src"),
    ):
        if os.path.isdir(os.path.join(cand, "models")):
            return cand
    return ""


def find_project_data() -> str:
    """定位本项目数据目录（含 conservation_profile.json），找不到返回 ''。"""
    for cand in (
        os.path.join(PROJECT_ROOT, "gv", "design", "data"),
        os.path.join(PROJECT_ROOT, "module1", "outputs"),
        os.path.join(PROJECT_ROOT, "module2", "out"),
    ):
        if os.path.isfile(os.path.join(cand, "conservation_profile.json")):
            return cand
    return ""


def ensure_mask_path():
    """返回本项目正式掩码路径（mod1 或 module4 official_mask.txt），存在则返回，否则返回 ''。"""
    for cand in (os.path.join(PROJECT_ROOT, "module1", "outputs", "mask.txt"),
                 os.path.join(PROJECT_ROOT, "gv", "design", "src", "module4", "config", "official_mask.txt")):
        if os.path.isfile(cand):
            return cand
    return ""


CANDIDATE_SRC = [
    os.environ.get("GVP_SRC", ""),
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"),
    find_project_src(),
    r"D:/aineeding/新建文件夹/gv/design/src",
    r"D:/aineeding/gv/design/src",
]


def ensure_file(path, desc):
    """确保文件存在，否则抛 FileNotFoundError 并给出明确提示。"""
    if not os.path.exists(path):
        raise FileNotFoundError("找不到%s: %s" % (desc, path))
    return path


def bootstrap_src() -> str:
    """把 gv 设计仓库 src/ 加入 sys.path，返回其绝对路径（找不到返回 ''）。

    优先取 GVP_SRC 环境变量；否则在若干候选路径中找含 models/ 的 src/。
    """
    for raw in CANDIDATE_SRC:
        base = os.path.expanduser(raw)
        models = os.path.join(base, "models")
        if base and os.path.isdir(models) and base not in sys.path:
            sys.path.insert(0, base)
            return base
    return ""


def load_optional_module(dotted, fromlist):
    """尽力导入 dotted 模块，成功返回模块对象，失败返回 None（不打断流程）。"""
    for mod in (dotted, dotted.rsplit(".", 1)[-1]):
        try:
            return __import__(mod, fromlist=fromlist)
        except Exception:
            continue
    return None


def read_fasta(path):
    """读取 fasta，返回 [(id, seq)]，序列保留原始大小写（按原文）。"""
    ensure_file(path, "输入文件")
    records, name, buf = [], None, []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name is not None:
        records.append((name, "".join(buf)))
    return records


def iter_sequence_lines(path):
    """逐行读取序列文件，返回去空白后的非空行列表。"""
    ensure_file(path, "文件")
    with open(path, "r", encoding="utf-8") as f:
        return [ln.strip() for ln in f if ln.strip()]


def read_locked_positions(path):
    """从掩码文件读 0 基锁定位点索引；支持 '#' 注释与整行索引。"""
    idx = []
    if not os.path.exists(path):
        return idx
    for ln in open(path, "r", encoding="utf-8"):
        s = ln.split("#")[0].strip()
        if not s:
            continue
        try:
            idx.append(int(s))
        except ValueError:
            pass
    return idx


def load_profile(data_dir, profile_name):
    """加载 conservation_profile.json。"""
    ensure_file(os.path.join(data_dir, profile_name), "conservation_profile.json")
    with open(os.path.join(data_dir, profile_name), "r", encoding="utf-8") as f:
        return json.load(f)


def write_csv(path, rows):
    """写 CSV，字段取所有行的并集，带 utf-8-sig（Excel 友好）。"""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    keys = sorted(set().union(*[set(r.keys()) for r in rows])) if rows else []
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def write_report(path, lines):
    """写入 UTF-8 文本报告。"""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
