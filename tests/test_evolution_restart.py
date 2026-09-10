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
        columns=reference['columns']
        for values in reference['history']:
            r=dict(zip(columns,values,strict=True))
            if abs(r['hydrogen_mass_Msun']-.1*r['central_X'])>2e-14:
                raise AssertionError('fully mixed hydrogen inventory differs from integrated fuel diagnostic')
            if not r['nuclear_deposited_Lsun']>0 or not r['nuclear_neutrino_Lsun']>0:
                raise AssertionError('nuclear luminosity diagnostics missing')
            if r['step_yr']==0:
                if r['last_halfstep_gravothermal_Lsun'] is not None:raise AssertionError('no half step exists at the initial state')
            elif abs((r['nuclear_deposited_Lsun']+r['last_halfstep_gravothermal_Lsun']-r['thermal_neutrino_Lsun'])/r['L_Lsun']-1-r['discrete_luminosity_balance'])>2e-14:
                raise AssertionError('reported energy components do not reproduce the accepted luminosity audit')
            if r['thermal_neutrino_Lsun']!=0:raise AssertionError('zero-loss control unexpectedly includes thermal neutrinos')
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
        # A long history must neither invalidate a checkpoint nor consume the
        # next invocation's step budget. Change counters in a test fixture only;
        # keep every physical variable and input identity exactly as written.
        lines=original_checkpoint.decode().splitlines()
        index=next(i for i,line in enumerate(lines) if line.startswith('128 '))
        header=lines[index].split();old_accepted=int(header[4]);old_rejected=int(header[5])
        header[4:6]=['15000','102'];lines[index]=' '.join(header)
        long_history=work/'long-history.restart'
        long_history.write_text('\n'.join(lines)+'\n')
        long_final=work/'long-final.restart'
        continued=run('long-history',['--restart',str(long_history),'--checkpoint',str(long_final)])
        if continued['profile']!=resumed['profile'] or continued['history']!=resumed['history']:
            raise AssertionError('lifetime counters changed the physical trajectory')
        if continued['rejected_steps']!=resumed['rejected_steps']-old_rejected+102:
            raise AssertionError('cumulative rejected steps were not preserved')
        roundtrip=run('long-final',['--restart',str(long_final)])
        if roundtrip['profile']!=reference['profile'] or roundtrip['restart']['accepted_steps_before_restart']!=15000+len(resumed['history'])-1:
            raise AssertionError('large lifetime counters did not round-trip')
        for invalid in ['-1','184467440737095516160']:
            header[4]=invalid;lines[index]=' '.join(header)
            bad=work/'invalid-counter.restart';bad.write_text('\n'.join(lines)+'\n')
            rejected=run('invalid-counter',['--restart',str(bad)],False)
            if 'counter' not in rejected['message']:raise AssertionError(rejected)
        rejected=run('changed-thermal-losses',['--restart',str(checkpoint),'--thermal-neutrinos','plasma-hrw'],False)
        if 'input tables differ' not in rejected['message']:raise AssertionError(rejected)
        plasma_checkpoint=work/'plasma.restart'
        plasma=run('plasma',['--thermal-neutrinos','plasma-hrw','--checkpoint',str(plasma_checkpoint),'--checkpoint-after','3'])
        if plasma['thermal_neutrino_model']!='plasma-hrw':raise AssertionError('missing thermal model selection')
        for values in plasma['history']:
            r=dict(zip(plasma['columns'],values,strict=True))
            if not r['thermal_neutrino_Lsun']>0:raise AssertionError('plasma cooling diagnostic missing')
            if r['step_yr'] and abs((r['nuclear_deposited_Lsun']+r['last_halfstep_gravothermal_Lsun']-r['thermal_neutrino_Lsun'])/r['L_Lsun']-1-r['discrete_luminosity_balance'])>2e-14:
                raise AssertionError('thermal sink missing from reported energy balance')
        plasma_resumed=run('plasma-resumed',['--thermal-neutrinos','plasma-hrw','--restart',str(plasma_checkpoint)])
        if plasma_resumed['profile']!=plasma['profile'] or plasma_resumed['history'][1:]!=plasma['history'][4:]:
            raise AssertionError('thermal-loss restart changed the accepted trajectory')
        explicit=run('explicit-original-opacity',['--opacity-directory',str(data/'opacity')])
        if explicit['history']!=reference['history'] or explicit['profile']!=reference['profile']:
            raise AssertionError('explicit original opacity changed the trajectory')
        # An alternate family with identical bytes is numerically equivalent,
        # but it is a different selection and must not reuse this checkpoint.
        alternate=work/'opacity';alternate.mkdir()
        for f in (data/'opacity').glob('*.dat'):shutil.copyfile(f,alternate/f.name)
        rejected=run('different-opacity-selection',
                     ['--restart',str(checkpoint),'--opacity-directory',str(alternate)],False)
        if 'input tables differ' not in rejected['message']:raise AssertionError(rejected)
        alternate_checkpoint=work/'alternate.restart'
        control=run('alternate-opacity',['--opacity-directory',str(alternate),
                                         '--checkpoint',str(alternate_checkpoint),'--checkpoint-after','3'])
        if control['history']!=reference['history'] or control['profile']!=reference['profile']:
            raise AssertionError('copied opacity family changed numerical results')
        # Verify that the selected plane bytes are pinned, not just the name
        # of the family. Whitespace leaves the parsed physics unchanged.
        plane=alternate/'tops_gs98_mixture_z020_low.dat'
        with plane.open('ab') as f:f.write(b'\n')
        rejected=run('changed-selected-opacity',
                     ['--restart',str(alternate_checkpoint),'--opacity-directory',str(alternate)],False)
        if 'input tables differ' not in rejected['message']:raise AssertionError(rejected)
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
                          'large_lifetime_counters_preserve_trajectory_and_roundtrip':True,
                          'rejects_changed_tolerance_changed_table_and_truncation':True}))


if __name__=='__main__':main()
