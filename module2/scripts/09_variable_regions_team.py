#!/usr/bin/env python3
"""Step 7：按团队口径输出「可变区域清单」（分工要求：预期约 37 位）。

分工文档要求模块2 交付「可变区域清单（预期约 37 位）」。该数字来自成员3 的
conservation_profile.json：majority_freq >= 0.9 锁定 35 个保守位 → 72 - 35 = 37 位可变。
本脚本把这 37 位逐位列出，并叠加本模块的熵/标签与 PF00741 域信息，供模块 3/5 直接使用。

输出：
    out/variable_regions_team_definition.csv   37 位清单（团队口径）
    out/variable_regions_summary.json          统计与两条口径的差异说明
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, 'out')


def main():
    m3 = json.load(open(os.path.join(OUT, 'member3_mask_reference.json')))
    ref_seq = m3['member3_ref_seq']
    locked = set(m3['mask_1based'])
    mf = m3.get('majority_freq')
    prof = {int(r['ref_pos']): r for r in csv.DictReader(open(os.path.join(OUT, 'ref_position_profile.csv')))}
    masks = json.load(open(os.path.join(OUT, 'mask_candidates_v2.json')))['masks']
    top40 = set(masks['top40']['positions'])
    dom = json.load(open(os.path.join(OUT, 'domain_stats.json')))['reference']['domain']
    dom_range = set(range(dom['ali_from'], dom['ali_to'] + 1))

    rows = []
    for p in range(1, len(ref_seq) + 1):
        if p in locked:
            continue
        pr = prof[p]
        rows.append({
            'ref_pos': p,
            'ref_residue': ref_seq[p - 1],
            'in_PF00741_domain': 'yes' if p in dom_range else 'no',
            'my_col_type': pr['col_type'],
            'my_entropy': pr['entropy'],
            'my_label': pr['label'],
            'member3_majority_freq': ('%.3f' % mf[p - 1]) if mf else '',
            'locked_by_my_top40': 'yes' if p in top40 else 'no',
            'note': ('域内高变位点，熵法同为 variable（两法一致）'
                     if pr['label'] == 'variable'
                     else '域外/侧翼，熵法为 insert（两法一致判为可变）'
                     if pr['label'] == 'insert'
                     else '域内中等变异：频率法判为可变，熵法为 intermediate（口径差异，建议 E2 敏感性验证）'),
        })

    with open(os.path.join(OUT, 'variable_regions_team_definition.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    summary = {
        'definition': '可变区域 = 参考序列 72 位中，未被成员3 掩码(majority_freq>=0.9, %d 位)锁定的位点'
                      % len(locked),
        'n_variable': len(rows),
        'n_in_domain': sum(1 for r in rows if r['in_PF00741_domain'] == 'yes'),
        'n_outside_domain': sum(1 for r in rows if r['in_PF00741_domain'] == 'no'),
        'n_also_locked_by_my_top40': sum(1 for r in rows if r['locked_by_my_top40'] == 'yes'),
        'consistency': {
            'consistent_variable(entropy also says variable or insert)': sum(
                1 for r in rows if r['my_label'] in ('variable', 'insert')),
            'definition_gap(entropy says intermediate/conserved)': sum(
                1 for r in rows if r['my_label'] in ('intermediate', 'conserved')),
        },
        'my_label_breakdown': {
            'variable': sum(1 for r in rows if r['my_label'] == 'variable'),
            'intermediate': sum(1 for r in rows if r['my_label'] == 'intermediate'),
            'conserved': sum(1 for r in rows if r['my_label'] == 'conserved'),
            'insert': sum(1 for r in rows if r['my_label'] == 'insert'),
        },
        'comparison_with_entropy_definition': {
            'this_file(team definition)': len(rows),
            'out/variable_positions.csv(entropy definition)': sum(
                1 for _ in open(os.path.join(OUT, 'variable_positions.csv'))) - 1,
            'note': '两条口径都保留：团队口径用于对齐成员3/模块3；熵口径用于 E2 掩码敏感性',
        },
        'source': 'out/member3_mask_reference.json + out/ref_position_profile.csv + out/mask_candidates_v2.json',
    }
    json.dump(summary, open(os.path.join(OUT, 'variable_regions_summary.json'), 'w'),
              ensure_ascii=False, indent=1)

    print('团队口径可变区域: %d 位（域内 %d / 域外 %d）' % (
        summary['n_variable'], summary['n_in_domain'], summary['n_outside_domain']))
    print('两法一致判为可变: %d 位；口径差异位点: %d 位' % (
        summary['consistency']['consistent_variable(entropy also says variable or insert)'],
        summary['consistency']['definition_gap(entropy says intermediate/conserved)']))
    print('-> out/variable_regions_team_definition.csv, out/variable_regions_summary.json')


if __name__ == '__main__':
    sys.exit(main())
