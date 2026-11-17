#!/usr/bin/env python3
"""Step 3：用 PF00741 对 856 条 GvpA 做 HMMER 域标注（hmmsearch, --cut_ga）。

输出：
    out/PF00741_hits.csv    每条序列的最佳域（ali/env 区间、E-value、score、覆盖率）
    out/domain_stats.json   汇总统计（命中率、域长分布、参考序列域区间）
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
FASTA = '/home/caizhi666/ml/pan_resources/gv/design/data/GvpA.fasta'
HMM = os.path.join(ROOT, 'pfam', 'PF00741.hmm')
OUT = os.path.join(ROOT, 'out')
REF_ID = os.environ.get('GVP_REF_ID', 'WP_007357596.1')  # 团队参考序列（与成员3 的 conservation_profile.json 一致）
CUTOFF = 'gathering'


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


def main():
    os.makedirs(OUT, exist_ok=True)
    recs = read_fasta(FASTA)
    seqlen = dict(recs)
    abc = Alphabet.amino()
    with HMMFile(HMM) as hf:
        hmm = hf.read()

    dig = [TextSequence(name=n.encode(), sequence=s.encode()).digitize(abc) for n, s in recs]
    # cpus=1 -> pyhmmer 走串行路径，避免多进程在沙箱内受限
    hits = list(pyhmmer.hmmsearch(hmm, dig, bit_cutoffs=CUTOFF, cpus=1))[0]

    rows = []
    for h in hits:
        if not len(h.domains):
            continue
        d = h.domains[0]
        L = len(seqlen.get(h.name, ''))
        rows.append({
            'seq_id': h.name,
            'seq_len': L,
            'ali_from': d.alignment.target_from,
            'ali_to': d.alignment.target_to,
            'env_from': d.env_from,
            'env_to': d.env_to,
            'evalue': '%.3g' % float(h.evalue),
            'score': round(float(h.score), 2),
            'dom_len': d.alignment.target_to - d.alignment.target_from + 1,
            'coverage': round((d.alignment.target_to - d.alignment.target_from + 1) / L, 4),
        })

    with open(os.path.join(OUT, 'PF00741_hits.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    spans = [r['dom_len'] for r in rows]
    covs = [r['coverage'] for r in rows]
    ref_rows = [r for r in rows if r['seq_id'] == REF_ID]
    stats = {
        'hmm': {'file': HMM, 'name': 'PF00741', 'M': hmm.M, 'cutoff': CUTOFF},
        'n_seq_total': len(recs),
        'n_seq_hit': len(rows),
        'hit_rate': round(len(rows) / len(recs), 4),
        'dom_len': {'min': min(spans), 'median': statistics.median(spans), 'max': max(spans)},
        'coverage': {'min': min(covs), 'median': round(statistics.median(covs), 4), 'max': max(covs)},
        'evalue_median': '%.3g' % statistics.median([float(r['evalue']) for r in rows]),
        'score_median': round(statistics.median([r['score'] for r in rows]), 2),
        'reference': {'id': REF_ID, 'len': len(seqlen.get(REF_ID, '')),
                      'domain': ref_rows[0] if ref_rows else None},
    }
    json.dump(stats, open(os.path.join(OUT, 'domain_stats.json'), 'w'), ensure_ascii=False, indent=1)

    print('PF00741 命中: %d/%d (cutoff=%s)' % (len(rows), len(recs), CUTOFF))
    print('域长: min %d 中位 %.0f max %d | 覆盖中位 %.2f' % (
        stats['dom_len']['min'], stats['dom_len']['median'], stats['dom_len']['max'],
        stats['coverage']['median']))
    if ref_rows:
        r = ref_rows[0]
        print('参考 %s: len=%d 域 %d-%d (E=%s, score=%s)' % (
            REF_ID, r['seq_len'], r['ali_from'], r['ali_to'], r['evalue'], r['score']))
    print('-> out/PF00741_hits.csv, out/domain_stats.json')


if __name__ == '__main__':
    sys.exit(main())
