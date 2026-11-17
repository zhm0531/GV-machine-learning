#!/usr/bin/env python3
"""Step 3b：用 PF00741 做 HMM 引导比对（hmmalign），作为保守度图谱的**主比对**。

为什么不用 FAMSA 做主线：
    scripts/02_build_msa.py 的 FAMSA 结果在 191 列中有 93 列是 100% gap（对短且高度
    分歧的 GvpA 家族过比对），不适合直接做逐列保守度。hmmalign 以 Pfam 家族 HMM 为
    骨架，输出 153 列且**无全 gap 列**，其中 39 列为 HMM 匹配列（高置信），其余为插入列。

输出：
    out/GvpA_hmmalign.sto        Stockholm 原始输出（保留大小写：大写=匹配列）
    out/GvpA_hmmalign.aln        同一比对的 FASTA 形式（保留大小写）
    out/GvpA_hmmalign_trim.aln   仅域内（trim=True）
    out/msa_hmmalign_meta.json   元信息
"""
import io
import json
import os
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


def write_fasta(path, names, rows, width=60):
    with open(path, 'w') as f:
        for n, s in zip(names, rows):
            f.write('>%s\n' % n)
            for i in range(0, len(s), width):
                f.write(s[i:i + width] + '\n')


def to_str(x):
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)


def main():
    os.makedirs(OUT, exist_ok=True)
    recs = read_fasta(FASTA)
    abc = Alphabet.amino()
    with HMMFile(HMM) as hf:
        hmm = hf.read()
    dig = [TextSequence(name=n.encode(), sequence=s.encode()).digitize(abc) for n, s in recs]

    msa = pyhmmer.hmmalign(hmm, dig, cpus=1, trim=False, all_consensus_cols=True)
    names = [to_str(n) for n in msa.names]
    rows = [to_str(a) for a in msa.alignment]
    ncol = len(rows[0])

    for p, label in ((os.path.join(OUT, 'GvpA_hmmalign.sto'), 'stockholm'),
                     (os.path.join(OUT, 'GvpA_hmmalign.a2m'), 'a2m')):
        buf = io.BytesIO()          # pyhmmer 的 writer 需要二进制流
        msa.write(buf, label)
        open(p, 'w').write(buf.getvalue().decode('utf-8', 'replace'))

    # 固定宽度 FASTA（保留大小写、含全部 gap 列），供下游按列解析
    write_fasta(os.path.join(OUT, 'GvpA_hmmalign_full.fa'), names, rows)

    msa_trim = pyhmmer.hmmalign(hmm, dig, cpus=1, trim=True, all_consensus_cols=True)
    write_fasta(os.path.join(OUT, 'GvpA_hmmalign_trim.aln'),
                [to_str(n) for n in msa_trim.names], [to_str(a) for a in msa_trim.alignment])

    gapfrac = [sum(1 for s in rows if s[i] in '.-') / len(rows) for i in range(ncol)]
    match_cols = [i for i in range(ncol)
                  if any(rows[j][i].isupper() for j in range(len(rows)))]
    meta = {
        'primary_alignment': 'hmmalign (pyhmmer) to PF00741',
        'tool': 'pyhmmer hmmalign',
        'n_seq': len(rows),
        'n_cols': ncol,
        'n_match_cols': len(match_cols),
        'n_allgap_cols': sum(1 for g in gapfrac if g > 0.999),
        'n_cols_gap_lt_0.5': sum(1 for g in gapfrac if g < 0.5),
        'trim_alignment_cols': len(to_str(msa_trim.alignment[0])),
        'reference': REF_ID,
        'cross_check_famsa': {
            'note': 'FAMSA 输出用于对照，未作为主比对',
            'file': os.path.join(OUT, 'GvpA_famsa.aln'),
        },
        'outputs': ['GvpA_hmmalign.sto', 'GvpA_hmmalign.a2m', 'GvpA_hmmalign_full.fa', 'GvpA_hmmalign_trim.aln'],
    }
    json.dump(meta, open(os.path.join(OUT, 'msa_hmmalign_meta.json'), 'w'),
              ensure_ascii=False, indent=1)

    print('主比对(hmmalign): %d 条 × %d 列（匹配列 %d，全gap列 %d）' % (
        len(rows), ncol, len(match_cols), meta['n_allgap_cols']))
    print('域内比对(trim): %d 列' % meta['trim_alignment_cols'])
    print('-> out/GvpA_hmmalign.{sto,a2m}, out/GvpA_hmmalign_trim.aln')


if __name__ == '__main__':
    sys.exit(main())
