#!/usr/bin/env python3
"""通用网盘下载器：从天大云盘分享链接按相对路径下载任意文件。

用法：
    python3 fetch_pan.py gv/recognition/data/picture/GvpA_region_16_46_five_panel.png ...
    python3 fetch_pan.py --out /some/dir path1 path2

原理：POST /link（取件码）→ link_token → /api/efast/v1/dir/list 逐层解析 docid/rev
      → /api/efast/v1/file/osdownload 取 S3 预签名 → GET 下载。
"""
import argparse
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = 'https://pan.tju.edu.cn'
LINK_ID = 'AAC85903DEA75D4DC3B8E604FED612552E'
PICKUP = 'oc7B'
ROOT_DOCID = ('gns://C30DE0CB64014B729A9BCF1AD0F20847/50A62FBC244C4B41BC7262FB7E9BAF7E/'
              '1AD1D810C9044AE197CBC2BA63735A8C/15FE6F023A8B4B679BC3C4A3E0CC8F2C')

TOKEN = None


def get_token():
    data = urllib.parse.urlencode({
        'id': LINK_ID, 'type': 'anonymous', 'password_required': 'true',
        'title': 'resources', 'item_type': 'folder', 'belongs_to': 'document',
        'password': PICKUP,
    }).encode()

    class NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, *a, **kw):
            return None

    opener = urllib.request.build_opener(NoRedirect)
    try:
        opener.open(urllib.request.Request(BASE + '/link', data=data), timeout=30)
    except urllib.error.HTTPError as e:
        sc = e.headers.get('Set-Cookie', '')
        m = re.search(r'link_token:' + LINK_ID + r'=([^;]+)', sc)
        if m:
            return m.group(1)
        raise RuntimeError('no token: ' + sc[:200])
    raise RuntimeError('unlock did not 302')


def api(path, payload):
    h = {'Authorization': 'Bearer ' + TOKEN, 'Accept': 'application/json',
         'Content-Type': 'application/json'}
    req = urllib.request.Request(BASE + path, data=json.dumps(payload).encode(),
                                 headers=h, method='POST')
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode('utf-8', 'replace'))


def children(docid):
    d = api('/api/efast/v1/dir/list', {'docid': docid, 'sort': 'asc', 'by': 'name'})
    out = {}
    for f in d.get('dirs', []):
        out[f['name']] = {'docid': f['docid'], 'is_dir': True}
    for f in d.get('files', []):
        out[f['name']] = {'docid': f['docid'], 'rev': f.get('rev'),
                          'size': f.get('size'), 'is_dir': False}
    return out


def resolve(relpath):
    cur = {'docid': ROOT_DOCID}
    for part in relpath.split('/'):
        kids = children(cur['docid'])
        if part not in kids:
            raise KeyError('未找到: ' + part)
        cur = kids[part]
    return cur


def download(entry, dest):
    r = api('/api/efast/v1/file/osdownload',
            {'docid': entry['docid'], 'rev': entry['rev'], 'offset': 0, 'length': -1})
    method, url, auth, amzdate = r['authrequest'][:4]
    hdrs = {}
    for h in (auth, amzdate):
        k, _, v = h.partition(':')
        hdrs[k.strip()] = v.strip()
    os.makedirs(os.path.dirname(dest) or '.', exist_ok=True)
    req = urllib.request.Request(url, headers=hdrs, method=method)
    with urllib.request.urlopen(req, timeout=300) as resp, open(dest, 'wb') as f:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
    return os.path.getsize(dest)


def main():
    global TOKEN
    ap = argparse.ArgumentParser()
    ap.add_argument('paths', nargs='+')
    ap.add_argument('--out', default='/home/caizhi666/ml/module2/reference')
    a = ap.parse_args()
    TOKEN = get_token()
    ok = fail = 0
    for rel in a.paths:
        dest = os.path.join(a.out, os.path.basename(rel))
        try:
            size = download(resolve(rel), dest)
            ok += 1
            print('OK   %-70s %8d B' % (rel, size))
        except Exception as e:
            fail += 1
            print('FAIL %-70s %s' % (rel, str(e)[:150]))
        time.sleep(0.2)
    print('完成: 成功 %d, 失败 %d' % (ok, fail))
    return 0 if fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
