#!/usr/bin/env bash
# 模块2 一键复现：从原始数据到全部交付物
#
# 用法：
#   cd module2 && bash scripts/run_all.sh
#
# 依赖（已在本目录 .venv 中装好，基于 uv）：
#   pyhmmer（HMMER3 的 Python 绑定）、pyfamsa、matplotlib
# 说明：
#   - pyhmmer 必须用 cpus=1 串行调用（沙箱下多进程共享内存不可用）
#   - 需要联网仅用于下载 Pfam HMM（PF00741 / PF01304 / PF05465），若 pfam/ 下已有可跳过
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="$ROOT/.venv/bin/python"
export TMPDIR="$ROOT/.tmp"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/home/caizhi666/ml/.mplcache}"
# 参考序列：默认与成员3 的 conservation_profile.json 一致（团队口径）
export GVP_REF_ID="${GVP_REF_ID:-WP_007357596.1}"
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$ROOT/out" "$ROOT/figures" "$ROOT/pfam"
echo "[i] 参考序列 GVP_REF_ID=$GVP_REF_ID"

# ---- 0) 下载 Pfam HMM（若缺失）----
fetch_hmm () {
  local acc="$1"
  if [ -f "$ROOT/pfam/$acc.hmm" ]; then return; fi
  echo "[0] 下载 $acc.hmm"
  curl -sS -L "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/$acc/?annotation=hmm" \
    -o "/tmp/$acc.hmm.gz"
  gunzip -c "/tmp/$acc.hmm.gz" > "$ROOT/pfam/$acc.hmm"
}
fetch_hmm PF00741      # GvpA 家族
fetch_hmm PF01304      # GvpC 重复单元（细菌）
fetch_hmm PF05465      # GvpC 重复单元（古菌，阴性对照）

# ---- 1) 域标注 ----
echo "[1] hmmsearch PF00741 -> out/PF00741_hits.csv"
"$PY" "$HERE/03_hmm_annotate.py"

# ---- 2) HMM 引导比对 ----
echo "[2] hmmalign -> out/GvpA_hmmalign_full.fa"
"$PY" "$HERE/04_hmmalign.py"

# ---- 2b) FAMSA 对照（可选，失败不影响主线）----
if [ ! -f "$ROOT/out/GvpA_famsa.aln" ]; then
  echo "[2b] FAMSA 对照比对（仅作对照，不作为主比对）"
  "$PY" "$HERE/02_build_msa.py" || echo "  FAMSA 失败，跳过（不影响主交付物）"
fi

# ---- 3) 保守度图谱 / 区间表 / 掩码候选 ----
echo "[3] 保守度图谱 -> out/conservation_profile.json 等"
"$PY" "$HERE/05_conservation_profile.py"

# ---- 4) GvpC 接口位点 ----
echo "[4] GvpC 重复标注与接口位点表 -> out/domain_interface_sites.csv"
"$PY" "$HERE/06_gvpc_interface.py"

# ---- 5) 图 ----
echo "[5] 生成 figures/*.png"
"$PY" "$HERE/07_figures.py"

# ---- 6) 与成员3 掩码对账（需要 out/member3_mask_reference.json，已随包提供）----
if [ -f "$ROOT/out/member3_mask_reference.json" ]; then
  echo "[6] 与 member3 掩码对账 -> out/reconcile_member3.json"
  "$PY" "$HERE/08_reconcile_member3.py"
  echo "[7] 团队口径可变区域清单（37 位）-> out/variable_regions_team_definition.csv"
  "$PY" "$HERE/09_variable_regions_team.py"
else
  echo "[6] 跳过对账与团队口径清单：缺少 out/member3_mask_reference.json"
  echo "    可从成员3 分支提取：git show origin/GV02-04-member3:gv/design/data/conservation_profile.json"
fi

echo
echo "全部完成。交付物见 out/ 与 figures/（清单见 README.md）"
