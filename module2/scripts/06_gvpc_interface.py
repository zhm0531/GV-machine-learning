#!/usr/bin/env python3
"""Step 4：GvpC 重复单元标注 + GvpA/GvpC 组装接口位点表。

方法：
  GvpC 侧 —— hmmsearch 比对 Pfam 重复域条目 PF01304（Gas_vesicle_C, LENG=33）
             与 PF05465（Halo_GVPC，古菌，作阴性对照），取 --cut_ga(gathering) 显著域。
  GvpA 侧 —— 取 PF00741 域区间作为与 GvpC 的结合面候选窗口。
  说明 —— 无实验复合物结构，残基级接触面属"推定"，全部标注置信度。

输出：
  out/gvpc_repeat_hits.csv        每条显著域（序列、区间、E-value、score）
  out/gvpc_repeat_summary.json    聚合统计（重复数、单元长度、覆盖比例）
  out/domain_interface_sites.csv  统一域/接口位点表（供模块3约束、模块5审计）
"""
import csv
import json
import os
import statistics
import sys

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import HMMFile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, 'out')
GVPC_FASTA = '/home/caizhi666/ml/pan_resources/gv/design/data/GvpC.fasta'
REF_ID = os.environ.get('GVP_REF_ID', 'WP_007357596.1')  # 团队参考序列（与成员3 的 conservation_profile.json 一致）
REPEAT_HMMS = ['PF01304', 'PF05465']


def read_fasta(path):
    recs, name, buf = [], None, []
    for line in open(path, encoding='utf-8', errors='replace'):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name is not None:
                recs.append((name, ''.join(buf)))
            name = line[1:].split()[0]
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        recs.append((name, ''.join(buf)))
    return recs


def hmm_domains(hmm_path, dig):
    with HMMFile(hmm_path) as hf:
        hmm = hf.read()
    hits = list(pyhmmer.hmmsearch(hmm, dig, bit_cutoffs='gathering', cpus=1))[0]
    doms = []
    for h in hits:
        for d in h.domains:
            doms.append({'seq_id': h.name, 'ali_from': d.alignment.target_from,
                         'ali_to': d.alignment.target_to,
                         'evalue': float(h.evalue), 'score': round(float(h.score), 2)})
    return doms, hmm


def main():
    recs = read_fasta(GVPC_FASTA)
    seqlen = {n: len(s) for n, s in recs}
    abc = Alphabet.amino()
    dig = [TextSequence(name=n.encode(), sequence=s.encode()).digitize(abc) for n, s in recs]

    all_doms, per_model = [], {}
    for acc in REPEAT_HMMS:
        doms, hmm = hmm_domains(os.path.join(ROOT, 'pfam', '%s.hmm' % acc), dig)
        per_model[acc] = {'n_domains': len(doms),
                          'n_seq_hit': len(set(d['seq_id'] for d in doms)),
                          'hmm_M': hmm.M}
        if acc == 'PF01304':
            all_doms = doms

    for d in all_doms:
        d['seq_len'] = seqlen.get(d['seq_id'], 0)
        d['dom_len'] = d['ali_to'] - d['ali_from'] + 1

    with open(os.path.join(OUT, 'gvpc_repeat_hits.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['seq_id', 'seq_len', 'ali_from', 'ali_to',
                                          'dom_len', 'evalue', 'score'])
        w.writeheader()
        w.writerows(all_doms)

    # 每条命中序列：重复块起点/终点/单元数
    by_seq = {}
    for d in all_doms:
        by_seq.setdefault(d['seq_id'], []).append(d)
    blocks = []
    for sid, doms in by_seq.items():
        doms.sort(key=lambda x: x['ali_from'])
        blocks.append({'seq_id': sid, 'seq_len': seqlen[sid], 'n_units': len(doms),
                       'block_start': doms[0]['ali_from'], 'block_end': doms[-1]['ali_to'],
                       'block_len': doms[-1]['ali_to'] - doms[0]['ali_from'] + 1,
                       'frac_of_seq': round((doms[-1]['ali_to'] - doms[0]['ali_from'] + 1)
                                            / seqlen[sid], 3),
                       'unit_len_median': statistics.median(
                           [d['dom_len'] for d in doms])})

    summary = {
        'gvpc': {'n_seq': len(recs),
                 'length_median': statistics.median(list(seqlen.values()))},
        'models': per_model,
        'n_seq': len(recs),
        'hit_sequences': len(blocks),
        'hit_rate': round(len(blocks) / len(recs), 4) if recs else 0,
        'units_per_seq_median': statistics.median([b['n_units'] for b in blocks]) if blocks else 0,
        'unit_len_median': statistics.median([b['unit_len_median'] for b in blocks]) if blocks else 0,
        'block_frac_of_seq_median': statistics.median([b['frac_of_seq'] for b in blocks]) if blocks else 0,
        'block_start_median': statistics.median([b['block_start'] for b in blocks]) if blocks else 0,
        'block_len_median': statistics.median([b['block_len'] for b in blocks]) if blocks else 0,
        'method': 'hmmsearch PF01304 --cut_ga(gathering); PF05465 作阴性对照',
    }
    json.dump(summary, open(os.path.join(OUT, 'gvpc_repeat_summary.json'), 'w'),
              ensure_ascii=False, indent=1)

    # ---- 统一域/接口位点表 ----
    ds = json.load(open(os.path.join(OUT, 'domain_stats.json')))
    ref_dom = ds['reference']['domain']
    sites = [
        {'site_id': 'GvpA_PF00741_domain', 'protein': 'GvpA(%s)' % REF_ID,
         'start': ref_dom['ali_from'], 'end': ref_dom['ali_to'], 'type': 'PF00741_domain',
         'evidence': 'hmmsearch PF00741 --cut_ga E=%s score=%s dom_len=%d (854/856 命中)'
                     % (ref_dom['evalue'], ref_dom['score'], ref_dom['dom_len']),
         'confidence': 'high', 'source': 'out/PF00741_hits.csv'},
        {'site_id': 'GvpA_GvpC_interface_candidate', 'protein': 'GvpA(%s)' % REF_ID,
         'start': ref_dom['ali_from'], 'end': ref_dom['ali_to'], 'type': 'assembly_interface',
         'evidence': '与 PF00741 域重合；GvpC 结合 GvpA 壳体表面为先验知识。'
                     '无实验复合物结构，残基级接触面未确定',
         'confidence': 'low', 'source': '推断（需结构/交联验证）'},
        {'site_id': 'GvpC_repeat_block', 'protein': 'GvpC(%d/%d 条有显著重复)'
                                                    % (summary['hit_sequences'], summary['n_seq']),
         'start': int(summary['block_start_median']),
         'end': int(summary['block_start_median'] + summary['block_len_median']),
         'type': 'assembly_interface',
         'evidence': 'PF01304 重复域：中位 %d 个单元 × %d aa，占序列 %.0f%%，块长中位 %d aa'
                     % (summary['units_per_seq_median'], summary['unit_len_median'],
                        100 * summary['block_frac_of_seq_median'], summary['block_len_median']),
         'confidence': 'medium', 'source': 'out/gvpc_repeat_hits.csv'},
        {'site_id': 'prior_work_regions', 'protein': 'GvpA/GvpC', 'start': 0, 'end': 0,
         'type': 'prior_work_reference',
         'evidence': '往届图：GvpA_region_16_46 / 25_55、GvpC_region_76_111 / 102_137'
                     '（module2/reference/，未复核，仅供参考）',
         'confidence': 'low', 'source': 'pan_resources reference figures'},
    ]
    with open(os.path.join(OUT, 'domain_interface_sites.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['site_id', 'protein', 'start', 'end', 'type',
                                          'evidence', 'confidence', 'source'])
        w.writeheader()
        w.writerows(sites)

    print('GvpC 重复标注(PF01304): 命中序列 %d/%d (%.1f%%)' % (
        summary['hit_sequences'], summary['n_seq'], 100 * summary['hit_rate']))
    print('  每序列单元数中位 %.0f | 单元长中位 %.0f aa | 重复块占序列 %.0f%% | 块长中位 %.0f aa' % (
        summary['units_per_seq_median'], summary['unit_len_median'],
        100 * summary['block_frac_of_seq_median'], summary['block_len_median']))
    print('  PF05465(古菌对照): 命中 %d' % per_model['PF05465']['n_seq_hit'])
    print('接口位点表: %d 行 -> out/domain_interface_sites.csv' % len(sites))


if __name__ == '__main__':
    sys.exit(main())
