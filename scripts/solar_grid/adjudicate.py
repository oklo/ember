"""Record a scoped scientific disposition and pin the evidence behind it."""
import argparse
import time
from pathlib import Path

from release import SPEC, digest, order_hash, orders, read, save


def record(root, key, judgement, reviewer, rationale, evidence):
    order = next((o for o in orders(root) if o['id'] == key), None)
    if order is None or not reviewer.strip() or not rationale.strip() or not evidence:
        raise ValueError('supply a known order, reviewer, rationale and evidence paths')
    pins = {}
    for name in evidence:
        source = (root / name).resolve()
        pins[str(source.relative_to(root.resolve()))] = digest(source)
    if judgement == 'adopted':
        states = read(root/'status.json')['recovery']
        if any(states.get(n) != 'passed' for n in order['targets']):
            raise ValueError('adoption requires evidence from every named target')
    value = dict(decision=judgement, reviewer=reviewer, rationale=rationale, evidence=pins,
                 at=time.time(), order_sha256=order_hash(order), release_sha256=digest(SPEC))
    save(root/'decisions'/key/(str(time.time_ns())+'.json'),value)
    save(root/'decisions'/(key+'.json'),value)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign',type=Path)
    parser.add_argument('order')
    parser.add_argument('decision',choices=['adopted','sensitivity','rejected'])
    parser.add_argument('--reviewer',required=True)
    parser.add_argument('--rationale',required=True)
    parser.add_argument('evidence',nargs='+')
    a=parser.parse_args()
    record(a.campaign, a.order, a.decision, a.reviewer, a.rationale, a.evidence)
