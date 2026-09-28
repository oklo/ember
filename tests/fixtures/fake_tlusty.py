#!/usr/bin/env python3
# Fake source for control-flow tests: reads stdin, writes complete fort.9 iterations from FAKE_HISTORY in ./fake.json
import json, sys, time, os
sys.stdin.read()
cfg = json.load(open('../../fake.json'))   # the column directory holds the scripted histories
orelax = 'ORELAX=1' in open('tas').read()
hist = cfg['newton' if orelax else 'damped']
with open('fort.9', 'w') as f:
    f.write(' RELATIVE CHANGES OF VECTOR PSI\n')
    for k, v in enumerate(hist, 1):
        for i in range(1, cfg['nd'] + 1):
            f.write(f'{k:5d}{i:5d}  {v:.2E}  0.00E+00  0.00E+00  0.00E+00  {v:.2E}    0    0\n')
        f.flush()
        with open('fort.7', 'w') as g:
            nd = cfg['nd']; g.write(f'{nd} -4\n' + '\n'.join(str(1e-6 * 10 ** (i / 4)) for i in range(nd)) + '\n'
                                    + '\n'.join('5000. 1e14 1e-8 1e17' for _ in range(nd)) + '\n')
        time.sleep(cfg.get('sleep', 0.6))
print('FAKE END')
