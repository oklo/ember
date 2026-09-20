#!/usr/bin/env python3
"""An interrupted numerical trajectory must match the uninterrupted trajectory.

Uses the actual stellar driver and tables, including non-equilibrium He3,
the metal EOS and non-grey boundary. No external Python packages are needed.
"""
import json
from pathlib import Path
import shutil
import shlex
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
        parallel_checkpoint=work/'parallel-midpoint.restart'
        parallel=run('parallel',['--step-workers','2','--checkpoint',str(parallel_checkpoint),
                                 '--checkpoint-after','3'])
        if parallel!=reference or parallel_checkpoint.read_bytes()!=checkpoint.read_bytes():
            raise AssertionError('two-worker evolution changed the trajectory or checkpoint bytes')
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
        parallel_resumed=run('parallel-resumed',['--restart',str(checkpoint),'--step-workers','2'],binary=copied)
        if parallel_resumed!=resumed:
            raise AssertionError('changing worker count across an exact restart changed the output')
        if checkpoint.read_bytes()!=original_checkpoint:raise AssertionError('restart input changed')
        if resumed['profile']!=reference['profile']:raise AssertionError('restart changed final stellar structure or composition')
        if resumed['history'][1:]!=reference['history'][4:]:raise AssertionError('restart changed the subsequent accepted trajectory')
        if resumed['restart']['history_start_age_yr']!=reference['history'][3][0]:raise AssertionError('restart age differs')
        if resumed['rejected_steps']!=reference['rejected_steps']:raise AssertionError('restart lost rejected-step count')
        # Synthetic extra rows are deliberately outside this short test track.
        # They test restart integrity, not the physical accuracy of new atmospheres.
        rows=original_table.decode().splitlines()
        ti=next(i for i,r in enumerate(rows) if r.startswith('log_teff '))
        gi=next(i for i,r in enumerate(rows) if r.startswith('log_g '))
        nt=int(rows[ti].split()[1]);ng=int(rows[gi].split()[1])
        top=float(rows[ti].split()[-1])
        rows[ti]='log_teff '+str(nt+1)+' '+' '.join(rows[ti].split()[2:])+f' {top+.05:.17g}'
        start=rows.index('data')+1
        expanded=[]
        for i in range(start,len(rows),nt*ng):
            plane=rows[i:i+nt*ng]
            if len(plane)!=nt*ng:raise AssertionError('unexpected fixture dimensions')
            expanded.extend(plane+plane[-ng:])
        warmer=work/'warmer-atmosphere.dat'
        warmer.write_text('\n'.join(rows[:start]+expanded)+'\n')
        warmer_arguments=arguments.copy();warmer_arguments[6]=f'nongrey:{warmer}'
        atmosphere_options=['--atmosphere-extension-restart',str(checkpoint),
                            '--restart-source-executable',str(executable),
                            '--restart-source-atmosphere',str(atmosphere)]
        warmer_checkpoint=work/'warmer.restart'
        warm=run('atmosphere-extension',atmosphere_options+['--checkpoint',str(warmer_checkpoint)],selection=warmer_arguments)
        if warm['history']!=resumed['history'] or warm['profile']!=resumed['profile']:
            raise AssertionError('appending unvisited atmosphere rows changed the resumed trajectory')
        if warm['atmosphere_extension']['added_states']!=len(expanded)//(nt+1):
            raise AssertionError('incorrect atmosphere extension state count')
        if checkpoint.read_bytes()!=original_checkpoint:
            raise AssertionError('atmosphere extension overwrote its source checkpoint')
        warm_roundtrip=run('warm-exact-restart',['--restart',str(warmer_checkpoint)],selection=warmer_arguments)
        if warm_roundtrip['profile']!=reference['profile']:
            raise AssertionError('new atmosphere identity did not survive exact restart')
        run('warm-without-extension',['--restart',str(checkpoint)],False,warmer_arguments)
        run('warm-wrong-source',atmosphere_options[:-1]+[str(warmer)],False,warmer_arguments)
        corrupt=warmer.read_text().replace('tau 100','tau 101')
        warmer.write_text(corrupt)
        rejected=run('warm-changed-depth',atmosphere_options,False,warmer_arguments)
        if 'matching depth' not in rejected['message']:raise AssertionError(rejected)
        warmer.write_text('\n'.join(rows[:start]+expanded)+'\n')
        # Extra EOS rows are synthetic and unvisited. Every original potential
        # coefficient is retained, including masked nodes and mixed derivatives.
        original_eos=data/'eos/freeeos300_gs98_z020.dat'
        hotter_eos=work/'hotter-eos';hotter_eos.mkdir()
        family_rows=original_eos.read_text().splitlines()
        hotter_family=hotter_eos/original_eos.name
        hotter_family.write_text('\n'.join(family_rows)+'\n')
        added_eos=0
        for name in family_rows[3:]:
            name=shlex.split(name)[0]
            contents=(original_eos.parent/name).read_text().splitlines()
            it=next(i for i,line in enumerate(contents) if line.startswith('log_t '))
            iq=next(i for i,line in enumerate(contents) if line.startswith('log_q '))
            axis=contents[it].split();nq=int(contents[iq].split()[1])
            contents[it]='log_t '+str(int(axis[1])+1)+' '+' '.join(axis[2:])+f' {float(axis[-1])+.025:.17g}'
            extra=contents[-nq:];added_eos+=sum(int(line.split()[0]) for line in extra)
            (hotter_eos/name).write_text('\n'.join(contents+extra)+'\n')
        hotter_arguments=arguments.copy();hotter_arguments[7]=f'metal:{hotter_family}'
        eos_options=['--eos-temperature-extension-restart',str(checkpoint),
                     '--restart-source-executable',str(executable),'--restart-source-eos',str(original_eos)]
        hotter_checkpoint=work/'hotter-eos.restart'
        hotter=run('eos-extension',eos_options+['--checkpoint',str(hotter_checkpoint)],selection=hotter_arguments)
        if hotter['history']!=resumed['history'] or hotter['profile']!=resumed['profile']:
            raise AssertionError('appending unvisited EOS rows changed the resumed trajectory')
        if hotter['eos_temperature_extension']['added_states']!=added_eos:
            raise AssertionError('incorrect added EOS state count')
        roundtrip=run('hotter-eos-exact',['--restart',str(hotter_checkpoint)],selection=hotter_arguments)
        if roundtrip['profile']!=reference['profile'] or checkpoint.read_bytes()!=original_checkpoint:
            raise AssertionError('EOS continuation did not preserve restart state')
        run('hotter-eos-without-extension',['--restart',str(checkpoint)],False,hotter_arguments)
        run('hotter-eos-wrong-source',eos_options[:-1]+[str(hotter_family)],False,hotter_arguments)
        changed_plane=hotter_eos/shlex.split(family_rows[3])[0]
        saved_plane=changed_plane.read_text();contents=saved_plane.splitlines();first=contents.index('data')+1
        fields=contents[first].split();fields[1]=format(float(fields[1])+1,'.17g');contents[first]=' '.join(fields)
        changed_plane.write_text('\n'.join(contents)+'\n')
        rejected=run('hotter-eos-changed-potential',eos_options,False,hotter_arguments)
        if 'potential values or masks changed' not in rejected['message']:raise AssertionError(rejected)
        contents=saved_plane.splitlines();fields=contents[first].split();fields[0]=str(1-int(fields[0]));contents[first]=' '.join(fields)
        changed_plane.write_text('\n'.join(contents)+'\n')
        rejected=run('hotter-eos-changed-mask',eos_options,False,hotter_arguments)
        if 'potential values or masks changed' not in rejected['message']:raise AssertionError(rejected)
        changed_plane.write_text(saved_plane.replace('composition_proxy "','composition_proxy "changed '))
        rejected=run('hotter-eos-changed-physics',eos_options,False,hotter_arguments)
        if 'source physics or composition changed' not in rejected['message']:raise AssertionError(rejected)
        changed_plane.write_text(saved_plane)
        # Density additions change row strides. Preserve every original node,
        # including the final column, and allow explicitly masked cold additions.
        denser_eos=work/'denser-eos';denser_eos.mkdir()
        denser_family=denser_eos/original_eos.name
        denser_family.write_text('\n'.join(family_rows)+'\n')
        added_density_states=0
        for name in family_rows[3:]:
            name=shlex.split(name)[0]
            contents=(original_eos.parent/name).read_text().splitlines()
            iq=next(i for i,line in enumerate(contents) if line.startswith('log_q '))
            axis=contents[iq].split();nq=int(axis[1]);start=contents.index('data')+1
            contents[iq]='log_q '+str(nq+1)+' '+' '.join(axis[2:])+f' {float(axis[-1])+.025:.17g}'
            expanded=[];nt=(len(contents)-start)//nq
            for it in range(nt):
                row=contents[start+it*nq:start+(it+1)*nq]
                extra=row[-1].split()
                if it<nt//2:extra[0]='0'
                added_density_states+=int(extra[0])
                expanded.extend(row+[' '.join(extra)])
            (denser_eos/name).write_text('\n'.join(contents[:start]+expanded)+'\n')
        denser_arguments=arguments.copy();denser_arguments[7]=f'metal:{denser_family}'
        density_options=['--eos-density-extension-restart',str(checkpoint),
                         '--restart-source-executable',str(executable),'--restart-source-eos',str(original_eos)]
        denser_checkpoint=work/'denser-eos.restart'
        denser=run('eos-density-extension',density_options+['--checkpoint',str(denser_checkpoint)],selection=denser_arguments)
        if denser['history']!=resumed['history'] or denser['profile']!=resumed['profile']:
            raise AssertionError('appending unvisited EOS columns changed the resumed trajectory')
        if denser['eos_density_extension']['added_states']!=added_density_states:
            raise AssertionError('incorrect valid density-extension state count')
        roundtrip=run('denser-eos-exact',['--restart',str(denser_checkpoint)],selection=denser_arguments)
        if roundtrip['profile']!=reference['profile'] or checkpoint.read_bytes()!=original_checkpoint:
            raise AssertionError('density continuation did not preserve restart state')
        run('denser-eos-without-extension',['--restart',str(checkpoint)],False,denser_arguments)
        run('denser-eos-temperature-option',eos_options,False,denser_arguments)
        changed_plane=denser_eos/shlex.split(family_rows[3])[0]
        saved_plane=changed_plane.read_text();contents=saved_plane.splitlines();start=contents.index('data')+1
        nq=int(next(line for line in contents if line.startswith('log_q ')).split()[1])
        index=start+2*nq-2  # Last original column of the second temperature row.
        fields=contents[index].split();fields[1]=format(float(fields[1])+1,'.17g');contents[index]=' '.join(fields)
        changed_plane.write_text('\n'.join(contents)+'\n')
        rejected=run('denser-eos-changed-original-column',density_options,False,denser_arguments)
        if 'potential values or masks changed' not in rejected['message']:raise AssertionError(rejected)
        changed_plane.write_text(saved_plane)
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
        # Synthetic lower-H planes exercise the extension contract only;
        # this short star stays inside the unchanged original X range.
        extended=work/'hot-extension';extended.mkdir()
        for f in (data/'opacity').glob('*.dat'):shutil.copyfile(f,extended/f.name)
        for f in extended.glob('tops_gs98_mixture_z*_high.dat'):
            rows=f.read_text().splitlines();header=rows[0].split();nx,nt,nr=map(int,header[:3])
            first=rows[3:4+nt];x,z=map(float,first[0].split())
            if x<=0:raise AssertionError('extension fixture requires a positive source floor')
            extra=first.copy();extra[0]=f'{x/2:.17g} {z:.17g}'
            header[0]=str(nx+1);rows[0]=' '.join(header)
            f.write_text('\n'.join(rows[:3]+extra+rows[3:])+'\n')
        extension_options=['--opacity-extension-restart',str(checkpoint),
                           '--restart-source-executable',str(executable),
                           '--restart-source-opacity',str(data/'opacity'),
                           '--opacity-directory',str(extended)]
        extended_checkpoint=work/'extended.restart'
        extension=run('hot-extension',extension_options+['--step-workers','2','--checkpoint',str(extended_checkpoint)])
        if extension['profile']!=reference['profile'] or extension['history']!=resumed['history']:
            raise AssertionError('opacity extension changed the original-domain trajectory')
        if extension['opacity_extension']['added_planes']!=3 or checkpoint.read_bytes()!=original_checkpoint:
            raise AssertionError('extension provenance missing or source checkpoint modified')
        extension_resumed=run('extension-exact-restart',['--restart',str(extended_checkpoint),'--opacity-directory',str(extended)])
        if extension_resumed['profile']!=reference['profile']:
            raise AssertionError('extended checkpoint does not resume exactly')
        f=extended/'tops_gs98_mixture_z020_high.dat';original_high=f.read_bytes()
        rows=original_high.decode().splitlines();nt=int(rows[0].split()[1]);j=5+nt
        values=rows[j].split();values[0]=format(float(values[0])+.01,'.17g');rows[j]=' '.join(values)
        f.write_text('\n'.join(rows)+'\n')
        rejected=run('changed-original-hot-entry',extension_options,False)
        if 'original opacity entries changed' not in rejected['message']:raise AssertionError(rejected)
        f.write_bytes(original_high)
        f=extended/'tops_gs98_mixture_z020_low.dat'
        with f.open('ab') as stream:stream.write(b'\n')
        rejected=run('changed-cool-extension',extension_options,False)
        if 'cooler opacity changed' not in rejected['message']:raise AssertionError(rejected)
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
