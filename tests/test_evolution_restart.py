#!/usr/bin/env python3
"""An interrupted numerical trajectory must match the uninterrupted trajectory.

Uses the actual stellar driver and tables, including non-equilibrium He3,
the metal EOS and non-grey boundary. No external Python packages are needed.
"""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    executable=Path(sys.argv[1]).resolve();data=Path(sys.argv[2]).resolve()
    with tempfile.TemporaryDirectory(prefix='ember-restart-test-') as temp:
        work=Path(temp);atmosphere=work/'atmosphere.dat'
        shutil.copyfile(data/'atmosphere/nongrey_gs98_z020_extended_tau100.dat',atmosphere)
        original_table=atmosphere.read_bytes()
        arguments=['128','1e9','1e7','10','sfii-svh','wd',f'nongrey:{atmosphere}',
                   f'metal:{data}/eos/freeeos300_gs98_z020.dat','ledoux-diffusive']
        def run(name,options,success=True,selection=None,binary=executable):
            result=subprocess.run([str(binary),*(selection or arguments),*options],capture_output=True,text=True)
            (work/(name+'.json')).write_text(result.stdout);(work/(name+'.log')).write_text(result.stderr)
            try:payload=json.loads(result.stdout)
            except ValueError:raise AssertionError((name,result.returncode,result.stdout,result.stderr))
            if (result.returncode==0)!=success or payload['converged']!=success:
                raise AssertionError((name,result.returncode,payload,result.stderr[-3000:]))
            return payload
        checkpoint=work/'midpoint.restart'
        reference=run('reference',['--checkpoint',str(checkpoint),'--checkpoint-after','3'])
        if not checkpoint.exists():raise AssertionError('requested checkpoint was not written')
        original_checkpoint=checkpoint.read_bytes()
        # Relocate the same executable, as the provenance runner does.
        copied=work/'ember-evolve';shutil.copy2(executable,copied)
        final=work/'final.restart'
        resumed=run('resumed',['--restart',str(checkpoint),'--checkpoint',str(final)],binary=copied)
        if checkpoint.read_bytes()!=original_checkpoint:raise AssertionError('restart input changed')
        if resumed['profile']!=reference['profile']:raise AssertionError('restart changed final stellar structure or composition')
        if resumed['history'][1:]!=reference['history'][4:]:raise AssertionError('restart changed the subsequent accepted trajectory')
        if resumed['restart']['history_start_age_yr']!=reference['history'][3][0]:raise AssertionError('restart age differs')
        if resumed['rejected_steps']!=reference['rejected_steps']:raise AssertionError('restart lost rejected-step count')
        finished=run('finished',['--restart',str(final)])
        if finished['profile']!=reference['profile']:raise AssertionError('final checkpoint did not round-trip exactly')
        changed=arguments.copy();changed[3]='20'
        rejected=run('different-tolerance',['--restart',str(checkpoint)],False,changed)
        if 'physics or tolerances differ' not in rejected['message']:raise AssertionError(rejected)
        atmosphere.write_bytes(original_table+b'\n')
        rejected=run('different-table',['--restart',str(checkpoint)],False)
        if 'input tables differ' not in rejected['message']:raise AssertionError(rejected)
        atmosphere.write_bytes(original_table)
        truncated=work/'truncated.restart';truncated.write_bytes(original_checkpoint[:len(original_checkpoint)//2])
        rejected=run('truncated',['--restart',str(truncated)],False)
        if 'checkpoint' not in rejected['message']:raise AssertionError(rejected)
        print(json.dumps({'accepted_trajectory_and_final_profile_identical':True,
                          'restart_age_yr':resumed['restart']['history_start_age_yr'],
                          'target_age_yr':reference['history'][-1][0],
                          'continued_steps':len(resumed['history'])-1,
                          'rejects_changed_tolerance_changed_table_and_truncation':True}))


if __name__=='__main__':main()
