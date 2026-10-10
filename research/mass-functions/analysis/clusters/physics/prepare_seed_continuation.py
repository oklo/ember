"""Private initializer experiment; evolution and stellar physics are unchanged.

Relax a compact model first, then increase the seed's imposed entropy loss in
small increments. These are static initialization solves, not elapsed evolution.
The final relaxed model alone starts the physical track at relative age zero.
This does not supply the missing early age or undo the post-deuterium assumption.
"""
from pathlib import Path
import hashlib, json, difflib

W = Path(__file__).resolve().parents[1]
source = W.parents[1] / 'out/luhman16-independent-oct8-v1/source/apps/lifetime_driver.cpp'
old = source.read_text()
new = old.replace(
    '    const double teff=cfg.number("initial_Teff_K"),entropy_loss=cfg.number("initial_entropy_loss");',
    '''    const double teff=cfg.number("initial_Teff_K"),entropy_loss=cfg.number("initial_entropy_loss");
    const double seed_start_loss=cfg.values.contains("initial_entropy_ramp_start")?
      cfg.number("initial_entropy_ramp_start"):entropy_loss;
    const double seed_ramp_steps=cfg.values.contains("initial_entropy_ramp_steps")?
      cfg.number("initial_entropy_ramp_steps"):0;
    if(!(seed_start_loss>0 && seed_ramp_steps>=0 && seed_ramp_steps<=100
         && seed_ramp_steps==std::floor(seed_ramp_steps)) ||
       (seed_ramp_steps==0 && seed_start_loss!=entropy_loss))
      throw std::invalid_argument("invalid initial entropy continuation");''')
assert new != old
new = new.replace('ContractingSource seed_source(nuclear,entropy_loss);',
                  'ContractingSource seed_source(nuclear,seed_start_loss);')
anchor = '      if(!solved.converged)throw std::runtime_error("Hayashi initial model: "+solved.message);'
assert new.count(anchor) == 1
new = new.replace(anchor, anchor + '''
      for(int stage=1;stage<=static_cast<int>(seed_ramp_steps);++stage) {
        seed_source.entropy_loss=seed_start_loss*std::pow(entropy_loss/seed_start_loss,
          static_cast<double>(stage)/seed_ramp_steps);
        auto next=relax(solved.model,seed_physics,atmosphere,options.relaxation);
        if(!next.converged)throw std::runtime_error("Initial entropy continuation stage "+
          std::to_string(stage)+": "+next.message);
        solved=std::move(next);
        std::cerr<<"Initial entropy continuation "<<stage<<" / "<<seed_ramp_steps
          <<" entropy_loss="<<seed_source.entropy_loss<<"\\n";
      }''')
dest = W / 'physics/lifetime_driver_seed_continuation.cpp'
dest.write_text(new)
(W/'physics/seed_continuation.diff').write_text(''.join(difflib.unified_diff(
    old.splitlines(True), new.splitlines(True), fromfile=str(source), tofile=str(dest))))
(W/'physics/seed_continuation_source.json').write_text(json.dumps(dict(
    base_commit='087351c', base_source=str(source),
    base_sha256=hashlib.sha256(old.encode()).hexdigest(),
    private_sha256=hashlib.sha256(new.encode()).hexdigest(), scope=__doc__,
), indent=2)+'\n')
