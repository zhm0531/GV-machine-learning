#!/usr/bin/env python3
"""Step 5：生成模块 2 的可视化图（3 张）。

  fig1_conservation_curve.png   逐列熵 + 滑动窗口均值 + 阈值 + 域/接口标注
  fig2_variability_heatmap.png  序列 × 列 的一致性热图（含 gap）
  fig3_region_map.png           保守/中间/可变/插入 区间图 + 域边界 + 掩码候选

图中文字使用英文，避免中文字体缺失导致方块。
"""
import csv
import json
import os
import random
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, 'out')
FIGS = os.path.join(ROOT, 'figures')

COLORS = {'conserved': '#2b6cb0', 'intermediate': '#90cdf4',
          'variable': '#e53e3e', 'insert': '#a0aec0'}


def read_full_aln(path):
    names, rows, name, buf = [], [], None, []
    for line in open(path, encoding='utf-8', errors='replace'):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name is not None:
                names.append(name)
                rows.append(''.join(buf))
            name = line[1:].split()[0]
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        names.append(name)
        rows.append(''.join(buf))
    return names, rows


def fig1(profile, dom, path):
    cols = profile['columns']
    x = [c['col'] for c in cols]
    y = [c['entropy'] for c in cols]
    w5 = profile['windows']['w5']
    q = profile['thresholds']['match_col_entropy_quantiles']

    fig, ax = plt.subplots(figsize=(15, 5.5))
    dom_cols = [c['col'] for c in cols if c['ref_pos'] and dom[0] <= c['ref_pos'] <= dom[1]]
    ax.axvspan(min(dom_cols), max(dom_cols), color='#fefcbf', alpha=0.8,
               label='PF00741 domain (ref %d-%d)' % dom)
    ax.bar(x, y, width=0.9, color='#cbd5e0', label='per-column entropy')
    for i, c in enumerate(cols):
        if c['label'] in COLORS:
            ax.bar([c['col']], [c['entropy']], width=0.9, color=COLORS[c['label']])
    ax.plot(x, w5, color='#1a202c', lw=1.8, label='sliding mean (w=5)')
    ax.axhline(q['25'], color='#2b6cb0', ls='--', lw=1, label='q25 (conserved cut)')
    ax.axhline(q['75'], color='#e53e3e', ls='--', lw=1, label='q75 (variable cut)')
    ax.set_xlabel('alignment column (hmmalign to PF00741)')
    ax.set_ylabel('Shannon entropy (bits)')
    ax.set_title('GvpA conservation profile (n=856; reference WP_073608715.1, 72 aa)')
    ax.legend(loc='upper right', fontsize=8, ncol=2)
    ax.set_xlim(0, len(cols) + 1)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig2(names, rows, profile, dom, path, n_show=60, seed=0):
    cols = profile['columns']
    cons = [c['consensus'] for c in cols]
    idx = list(range(len(rows)))
    random.Random(seed).shuffle(idx)
    idx = sorted(idx[:n_show])
    mat = []
    for i in idx:
        s = rows[i]
        mat.append([0.5 if s[j] in '.-' else (1.0 if s[j].upper() == cons[j] else 0.0)
                    for j in range(len(cols))])
    fig, ax = plt.subplots(figsize=(15, 6))
    ax.imshow(mat, aspect='auto', cmap='RdYlGn', vmin=0, vmax=1, interpolation='nearest')
    dom_cols = [c['col'] - 1 for c in cols if c['ref_pos'] and dom[0] <= c['ref_pos'] <= dom[1]]
    ax.axvline(min(dom_cols), color='black', lw=1.2)
    ax.axvline(max(dom_cols), color='black', lw=1.2)
    ax.set_xlabel('alignment column')
    ax.set_ylabel('sampled sequences (%d of %d)' % (n_show, len(rows)))
    ax.set_title('Residue identity to column consensus '
                 '(green=match, red=mismatch, yellow=gap; black lines = PF00741 domain)')
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig3(regions, profile, dom, path):
    fig, ax = plt.subplots(figsize=(15, 2.6))
    ncol = len(profile['columns'])
    for r in regions:
        ax.add_patch(Rectangle((r['start_col'], 0), r['end_col'] - r['start_col'] + 1, 1,
                               color=COLORS.get(r['label'], '#cbd5e0')))
    dom_cols = [c['col'] for c in profile['columns']
                if c['ref_pos'] and dom[0] <= c['ref_pos'] <= dom[1]]
    ax.axvline(min(dom_cols), color='black', lw=2, ls='--')
    ax.axvline(max(dom_cols), color='black', lw=2, ls='--')
    ax.text(min(dom_cols), 1.08, 'PF00741 domain', fontsize=9)
    ax.set_xlim(0, ncol)
    ax.set_ylim(0, 1.25)
    ax.set_yticks([])
    ax.set_xlabel('alignment column')
    ax.set_title('Conserved / intermediate / variable / insert regions')
    handles = [Rectangle((0, 0), 1, 1, color=COLORS[k]) for k in
               ['conserved', 'intermediate', 'variable', 'insert']]
    ax.legend(handles, ['conserved', 'intermediate', 'variable', 'insert (flank)'],
              ncol=4, fontsize=8, loc='lower right')
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    os.makedirs(FIGS, exist_ok=True)
    profile = json.load(open(os.path.join(OUT, 'conservation_profile.json')))
    ds = json.load(open(os.path.join(OUT, 'domain_stats.json')))
    dom = (ds['reference']['domain']['ali_from'], ds['reference']['domain']['ali_to'])
    regions = list(csv.DictReader(open(os.path.join(OUT, 'conservation_regions.csv'))))
    for r in regions:
        r['start_col'] = int(r['start_col'])
        r['end_col'] = int(r['end_col'])

    fig1(profile, dom, os.path.join(FIGS, 'fig1_conservation_curve.png'))
    names, rows = read_full_aln(os.path.join(OUT, 'GvpA_hmmalign_full.fa'))
    fig2(names, rows, profile, dom, os.path.join(FIGS, 'fig2_variability_heatmap.png'))
    fig3(regions, profile, dom, os.path.join(FIGS, 'fig3_region_map.png'))
    for f in sorted(os.listdir(FIGS)):
        print('-> figures/%s (%d B)' % (f, os.path.getsize(os.path.join(FIGS, f))))


if __name__ == '__main__':
    sys.exit(main())
