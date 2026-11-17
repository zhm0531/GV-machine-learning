#!/usr/bin/env python3
"""Step 0：用 FAMSA 对 856 条 GvpA 全长序列重建正式多序列比对（替代往届的截断伪比对）。

输出：
    out/GvpA_famsa.aln       对齐后的 FASTA
    out/msa_meta.json        工具版本 / 序列数 / 列数 / 耗时 / 输入指纹
"""
import hashlib
import json
import os
import sys
import time

from pyfamsa import Aligner, Sequence

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FASTA = '/home/caizhi666/ml/pan_resources/gv/design/data/GvpA.fasta'
OUT = os.path.join(ROOT, 'out')
THREADS = int(os.environ.get('THREADS', '8'))


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


def write_fasta(path, rows, width=60):
    with open(path, 'w') as f:
        for n, s in rows:
            f.write('>%s\n' % n)
            for i in range(0, len(s), width):
                f.write(s[i:i + width] + '\n')


def main():
    os.makedirs(OUT, exist_ok=True)
    recs = read_fasta(FASTA)
    if not recs:
        sys.exit('序列表为空: ' + FASTA)
    print('输入序列: %d 条' % len(recs))

    seqs = [Sequence(n.encode(), s.encode()) for n, s in recs]
    t0 = time.time()
    aligner = Aligner(threads=THREADS)
    aln = aligner.align(seqs)
    rows = [(g.id.decode(), g.sequence.decode().upper()) for g in aln]
    dt = time.time() - t0

    ncol = len(rows[0][1])
    assert all(len(s) == ncol for _, s in rows), '比对长度不一致'
    out_aln = os.path.join(OUT, 'GvpA_famsa.aln')
    write_fasta(out_aln, rows)

    meta = {
        'input_fasta': FASTA,
        'input_sha1': hashlib.sha1(open(FASTA, 'rb').read()).hexdigest()[:16],
        'n_seq': len(rows),
        'n_cols': ncol,
        'tool': 'pyfamsa (FAMSA)',
        'threads': THREADS,
        'elapsed_sec': round(dt, 2),
        'output': out_aln,
    }
    json.dump(meta, open(os.path.join(OUT, 'msa_meta.json'), 'w'), indent=1)
    print('比对完成: %d 条 × %d 列, 用时 %.1fs' % (len(rows), ncol, dt))
    print('-> %s' % out_aln)


if __name__ == '__main__':
    main()
