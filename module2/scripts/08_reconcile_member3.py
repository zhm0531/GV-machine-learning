#!/usr/bin/env python3
"""Step 6（对账）：把本模块的保守度结论与成员3 的 conservation_profile.json 逐位比对。

背景：分工文档要求模块2「复用 conservation_profile.json」——该文件是成员3 在 M2 的产出
（仓库分支 GV02-04-member3 的 gv/design/data/conservation_profile.json），
参考序列 WP_007357596.1，掩码口径 majority_freq >= 0.9 → 35 个锁定位。
本脚本把两条独立路径的结果对齐到同一参考坐标后比较，输出差异清单，避免下游（模块1/3）
拿到两套互相矛盾的掩码。

用法：
    python3 scripts/08_reconcile_member3.py                 # 使用 out/member3_mask_reference.json
    python3 scripts/08_reconcile_member3.py --member3 /path/to/conservation_profile.json

输出：
    out/reconcile_member3.json    对账指标
    figures/fig4_mask_comparison.png  两套掩码对比条带图
"""
import argparse
import csv
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, 'out')
FIGS = os.path.join(ROOT, 'figures')


def load_member3(path):
    """读取成员3 的 profile（或本包提取的精简版）。掩码统一转成 1-based。"""
    d = json.load(open(path))
    if 'mask_1based' in d:                      # 精简版
        return {'ref_id': d['member3_ref_id'], 'ref_seq': d['member3_ref_seq'],
                'mask': set(d['mask_1based']), 'thr': d.get('mask_threshold'),
                'majority_freq': d.get('majority_freq')}
    mask = d['mask']
    if mask and isinstance(mask[0], bool):
        pos = {i + 1 for i, v in enumerate(mask) if v}
    else:
        pos = {int(x) + 1 for x in mask}        # 原文件为 0-based
    return {'ref_id': d['ref_id'], 'ref_seq': d.get('ref_seq', ''),
            'mask': pos, 'thr': d.get('mask_threshold'),
            'majority_freq': d.get('majority_freq')}


def spearman(a, b):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    cov = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    va = sum((x - ma) ** 2 for x in ra) ** 0.5
    vb = sum((x - mb) ** 2 for x in rb) ** 0.5
    return cov / (va * vb) if va and vb else 0.0


def jaccard(a, b):
    return len(a & b) / len(a | b) if (a | b) else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--member3', default=os.path.join(OUT, 'member3_mask_reference.json'))
    a = ap.parse_args()

    m3 = load_member3(a.member3)
    prof = list(csv.DictReader(open(os.path.join(OUT, 'ref_position_profile.csv'))))
    for r in prof:
        r['ref_pos'] = int(r['ref_pos'])
        r['entropy'] = float(r['entropy'])
    by_pos = {r['ref_pos']: r for r in prof}
    masks = json.load(open(os.path.join(OUT, 'mask_candidates_v2.json')))['masks']
    dom = json.load(open(os.path.join(OUT, 'domain_stats.json')))['reference']['domain']
    dom_range = set(range(dom['ali_from'], dom['ali_to'] + 1))

    mine_conserved = {r['ref_pos'] for r in prof if r['label'] == 'conserved'}
    mine_variable = {r['ref_pos'] for r in prof if r['label'] == 'variable'}
    mine_entropy = [by_pos[p]['entropy'] for p in range(1, len(prof) + 1)]

    res = {
        'reference': {'member3_ref_id': m3['ref_id'], 'len': len(prof),
                      'member3_mask_threshold': m3['thr'],
                      'member3_n_locked': len(m3['mask']),
                      'member3_n_free': len(prof) - len(m3['mask'])},
        'member3_mask_1based': sorted(m3['mask']),
        'mine': {
            'n_conserved_label': len(mine_conserved),
            'n_variable_label': len(mine_variable),
            'domain_range': [dom['ali_from'], dom['ali_to']],
        },
        'overlap': {},
        'distribution': {
            'member3_locked_inside_domain': len(m3['mask'] & dom_range),
            'member3_locked_outside_domain': len(m3['mask'] - dom_range),
            'member3_locked_in_my_insert_cols': sum(
                1 for p in m3['mask'] if by_pos[p]['col_type'] == 'insert'),
        },
    }
    for name, mk in masks.items():
        s = set(mk['positions'])
        inter = m3['mask'] & s
        res['overlap'][name] = {
            'n_mine': len(s), 'n_member3': len(m3['mask']), 'n_intersect': len(inter),
            'jaccard': round(jaccard(m3['mask'], s), 4),
            'precision_vs_member3': round(len(inter) / len(s), 4) if s else 0,
            'recall_vs_member3': round(len(inter) / len(m3['mask']), 4) if m3['mask'] else 0,
        }
    inter = m3['mask'] & mine_conserved
    res['overlap']['my_conserved_label(q25)'] = {
        'n_mine': len(mine_conserved), 'n_member3': len(m3['mask']), 'n_intersect': len(inter),
        'jaccard': round(jaccard(m3['mask'], mine_conserved), 4),
        'precision_vs_member3': round(len(inter) / len(mine_conserved), 4) if mine_conserved else 0,
        'recall_vs_member3': round(len(inter) / len(m3['mask']), 4) if m3['mask'] else 0,
    }
    if m3['majority_freq'] and len(m3['majority_freq']) == len(prof):
        # 成员3 的 majority_freq 越高越保守；我的 entropy 越低越保守 → 期望强负相关
        res['spearman_majorityfreq_vs_entropy'] = round(
            spearman(m3['majority_freq'], [-e for e in mine_entropy]), 4)
    json.dump(res, open(os.path.join(OUT, 'reconcile_member3.json'), 'w'),
              ensure_ascii=False, indent=1)

    # ---- 对比图 ----
    n = len(prof)
    fig, ax = plt.subplots(figsize=(14, 2.6))
    for p in sorted(m3['mask']):
        ax.add_patch(Rectangle((p - 0.5, 1.05), 1, 0.8, color='#2b6cb0'))
    top = set(masks['top40']['positions'])
    for p in sorted(top):
        ax.add_patch(Rectangle((p - 0.5, 0.05), 1, 0.8, color='#e53e3e'))
    ax.axvspan(dom['ali_from'] - 0.5, dom['ali_to'] + 0.5, color='#fefcbf', zorder=0)
    ax.set_xlim(0.5, n + 0.5)
    ax.set_ylim(0, 2)
    ax.set_yticks([1.45, 0.45])
    ax.set_yticklabels(['member3 mask\n(majority_freq>=%.1f, n=%d)'
                        % (m3['thr'] or 0, len(m3['mask'])),
                        'this module top40\n(entropy rank, n=%d)' % len(top)], fontsize=8)
    ax.set_xlabel('reference position (%s, 1-based)' % m3['ref_id'])
    ax.set_title('Mask reconciliation: member3 vs this module (shaded = PF00741 domain %d-%d)'
                 % (dom['ali_from'], dom['ali_to']))
    fig.tight_layout()
    os.makedirs(FIGS, exist_ok=True)
    fig.savefig(os.path.join(FIGS, 'fig4_mask_comparison.png'), dpi=150)
    plt.close(fig)

    o = res['overlap']['top40']
    print('成员3 掩码: %d 锁定 / %d 自由（参考 %s）' % (
        len(m3['mask']), len(prof) - len(m3['mask']), m3['ref_id']))
    print('本模块 top40: %d 位，与成员3 交集 %d，Jaccard %.3f' % (
        o['n_mine'], o['n_intersect'], o['jaccard']))
    print('成员3 锁定位落在 PF00741 域内: %d，域外: %d' % (
        res['distribution']['member3_locked_inside_domain'],
        res['distribution']['member3_locked_outside_domain']))
    if 'spearman_majorityfreq_vs_entropy' in res:
        print('成员3 majority_freq 与我的熵（取负）Spearman: %.3f' % res['spearman_majorityfreq_vs_entropy'])
    print('-> out/reconcile_member3.json, figures/fig4_mask_comparison.png')


if __name__ == '__main__':
    sys.exit(main())
