#!/usr/bin/env python3
"""Source acceptance tests; synthetic records only test the import contract."""
import math
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import struct
import shlex
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from generate_nongrey_grid import composition, resample_initial_structure, input_fingerprint, continuation_structure, truncate_initial_structure, scale_initial_column, run_opacity_jobs, atmosphere_inputs, validate_source_capacity
from import_nongrey_grid import source_state, source_inputs, import_grid, source_diagnostics_match
from nongrey_opacity import read_table, validate_table, merge_isotherms
from assemble_nongrey_grid import complete_cells


class SourceAcceptance(unittest.TestCase):
    def assert_patch_digest(self, path, expected):
        # Archived calculations retain their own source version when the
        # current patch advances. Verify those bytes rather than relabel them.
        contents = path.read_bytes()
        if hashlib.sha256(contents).hexdigest() != expected:
            archive = path.parent/"nongrey_validation"/f"{expected}.patch.gz"
            contents = gzip.decompress(archive.read_bytes())
        self.assertEqual(hashlib.sha256(contents).hexdigest(), expected)

    def test_numerical_solver_comparison_requires_unchanged_evidence(self):
        from extend_mixed_atmosphere import numerical_variants
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory); source = path/'source'; source.write_text('reviewed source')
            report = dict(passed=True, maximum_allowed_relative_difference=.001,
                          variant_executable_sha256='new', reference_executable_sha256='old',
                          cases=[dict(accepted=True, relative_changes=dict(T=1e-7, Pgas=-2e-7))],
                          input_sha256={'source':hashlib.sha256(source.read_bytes()).hexdigest()})
            receipt = path/'comparison.json'; receipt.write_text(json.dumps(report))
            variants, evidence = numerical_variants([receipt])
            self.assertEqual(variants, {'new':'old'}); self.assertIn(str(receipt), evidence)
            source.write_text('changed')
            with self.assertRaisesRegex(ValueError, 'dependency changed'):
                numerical_variants([receipt])
            source.write_text('reviewed source'); report['cases'][0]['relative_changes']['Pgas']=.01
            receipt.write_text(json.dumps(report))
            with self.assertRaisesRegex(ValueError, 'comparison failed'):
                numerical_variants([receipt])

    def test_mixed_extension_preserves_values_and_disabled_cells(self):
        from extend_mixed_atmosphere import render
        axes = [[.98, .9955], [3600., 3800.], [6.3, 6.5]]
        header = ['EMBER_COMPOSITION_ATMOSPHERE 4', 'hydrogen 2 .98 .9955',
                  'helium3 1 0', 'log_teff 2 3.5563025007672873 3.5797835966168101',
                  'log_g 2 6.3 6.5', 'data']
        old = {key: f'1 {3.7+i*.001:.17g} {9+i*.01:.17g}'
               for i, key in enumerate(itertools.product(*axes))}
        cell = ((.98, .9955), (3600., 3800.), (6.3, 6.5))
        additions = {(h, 3400., g): dict(T=4800., Pgas=1e9)
                     for h, g in itertools.product(axes[0], axes[2])}
        text, new_axes, cells = render(header, axes, old, {cell: 0}, additions)
        lines = text.splitlines()
        rows = dict(zip(itertools.product(*new_axes),
                        lines[lines.index('data')+1:lines.index('cells')], strict=True))
        self.assertTrue(all(rows[k] == value for k, value in old.items()))
        self.assertEqual([c['supported'] for c in cells], [True, False])
        self.assertIn('hydrogen 2 .98 .9955', lines)
        self.assertEqual(new_axes[1], [3400., 3600., 3800.])

    def test_mixed_extension_rejects_changed_old_interpolation_and_partial_cells(self):
        from extend_mixed_atmosphere import render
        axes = [[.98, .9955], [3600., 3800.], [6.3, 6.5]]
        header = ['hydrogen 2 .98 .9955', 'log_teff 2 3.5563 3.5798',
                  'log_g 2 6.3 6.5', 'data']
        old = {k: '1 3.7 9' for k in itertools.product(*axes)}
        cell = ((.98, .9955), (3600., 3800.), (6.3, 6.5))
        state = dict(T=4800., Pgas=1e9)
        with self.assertRaisesRegex(ValueError, 'interior coordinate'):
            render(header, axes, old, {cell: 1}, {(.98, 3700., 6.3): state})
        with self.assertRaisesRegex(ValueError, 'replace an existing'):
            render(header, axes, old, {cell: 1}, {(.98, 3600., 6.3): state})
        _, _, cells = render(header, axes, old, {cell: 1}, {(.98, 3400., 6.3): state})
        self.assertEqual([c['supported'] for c in cells], [False, True])

    def test_mixed_extension_rejects_unfinished_source(self):
        from extend_mixed_atmosphere import load_column
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path/'result.json').write_text(json.dumps(dict(accepted=False)))
            with self.assertRaisesRegex(ValueError, 'unaccepted atmosphere'):
                load_column(path, {})

    def test_depth_control_requires_verified_donor_and_deeper_boundary(self):
        from generate_nongrey_grid import main
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prepared, spec = root/'prepared.json', root/'spec.json'
            prepared.write_text('{}')
            spec.write_text(json.dumps({'tau':100., 'depths':100}))
            base = ['generate_nongrey_grid.py', str(prepared), str(spec), str(root/'work')]
            invalid = [
                ['--initial-bottom-tau','300'],
                ['--initial-bottom-tau','100','--composition-donor',str(root/'donor')],
                ['--initial-bottom-tau','nan','--continuation-models',str(root/'donor')],
                ['--initial-bottom-tau','300','--composition-donor',str(root/'donor'),
                 '--initial-models',str(root/'other')],
            ]
            for arguments in invalid:
                with self.subTest(arguments=arguments), patch.object(sys,'argv',base+arguments):
                    with self.assertRaisesRegex(ValueError,'initial bottom depth'):
                        main()
            self.assertFalse((root/'work').exists())

    def test_capacity_requires_matching_build_and_enforces_all_dimensions(self):
        spec = dict(log_temperature=[31, 3.2, 4.2], log_density=[19, -13, -3],
                    depths=100, atmosphere_frequencies=20000, opacity_frequencies=30000)
        with self.assertRaisesRegex(ValueError, "compiled source capacity"):
            validate_source_capacity({}, spec)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binary, header, receipt = [root/p for p in ["tlusty", "BASICS.FOR", "build.json"]]
            binary.write_bytes(b"synthetic executable")
            header.write_text("C     PARAMETER (MTABT=21)\n"
                              "      PARAMETER (MTABT=35, MTABR=37,\n"
                              "     * MDEPTH=400, MFREQ=32000, MFRTAB=30005) ! selected build\n")
            sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
            receipt.write_text(json.dumps(dict(executable_sha256=sha(binary),
                                               input_sha256={str(header): sha(header)})))
            prepared = dict(tlusty=str(binary), executables={"tlusty": sha(binary)},
                            tlusty_build=str(receipt), tlusty_build_sha256=sha(receipt))
            self.assertEqual(validate_source_capacity(prepared, spec)["MTABT"], 35)
            for key, value in [("log_temperature", [36, 3.2, 4.2]),
                               ("log_density", [20, -13, -3]), ("depths", 401),
                               ("atmosphere_frequencies", 32001), ("opacity_frequencies", 30006)]:
                with self.subTest(key=key), self.assertRaises(ValueError):
                    validate_source_capacity(prepared, dict(spec, **{key: value}))
            binary.write_bytes(b"different executable")
            with self.assertRaisesRegex(ValueError, "selected executable"):
                validate_source_capacity(prepared, spec)
            binary.write_bytes(b"synthetic executable")
            header.write_text(header.read_text()+"C changed\n")
            with self.assertRaisesRegex(ValueError, "header missing or changed"):
                validate_source_capacity(prepared, spec)
            receipt.write_text(receipt.read_text()+"\n")
            with self.assertRaisesRegex(ValueError, "build record changed"):
                validate_source_capacity(prepared, spec)

    def test_explicit_relaxation_changes_only_numerical_input_and_fingerprint(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root/"data").mkdir()
            table = root/"table.bin"
            table.write_bytes(b"synthetic table")
            spec = dict(alpha=1.9, atmosphere_frequencies=5000, depths=100,
                        tau_top=1e-7, tau_bottom=1000.)
            prepared = dict(data_sha256="synthetic", synple=str(root))
            abundance, masses = composition(.999, 0, [0]*5)
            plain, damped = root/"plain", root/"damped"
            atmosphere_inputs(plain, prepared, spec, table, abundance, masses, 4600, 5.65)
            atmosphere_inputs(damped, prepared, dict(spec, newton_relaxation=.3),
                              table, abundance, masses, 4600, 5.65)
            self.assertNotIn("ORELAX", (plain/"tas").read_text())
            self.assertEqual((damped/"tas").read_text(),
                             (plain/"tas").read_text()+"ORELAX=0.3\n")
            for name in ["fort.5", "ember-masses.dat", "opacity.sha256"]:
                self.assertEqual((plain/name).read_bytes(), (damped/name).read_bytes())
            self.assertNotEqual(input_fingerprint("same-binary", plain),
                                input_fingerprint("same-binary", damped))
            tightened = root/"tightened"
            atmosphere_inputs(tightened, prepared, dict(spec, temperature_convergence=1e-8),
                              table, abundance, masses, 4600, 5.65)
            self.assertEqual((tightened/"tas").read_text(),
                             (plain/"tas").read_text().replace("CHMAX=1.e-6", "CHMAX=1e-08"))
            for name in ["fort.5", "ember-masses.dat", "opacity.sha256"]:
                self.assertEqual((plain/name).read_bytes(), (tightened/name).read_bytes())
            self.assertNotEqual(input_fingerprint("same-binary", plain),
                                input_fingerprint("same-binary", tightened))
            # .17g formerly expanded 1e-7 to a 22-character token, aborting
            # TLUSTY before it could read the atmosphere or solve any physics.
            for threshold in [1e-7, 1.2e-7]:
                changed = root/str(threshold)
                atmosphere_inputs(changed, prepared, dict(spec, temperature_convergence=threshold),
                                  table, abundance, masses, 4600, 5.65)
                token = (changed/"tas").read_text().split("CHMAX=")[1].replace(",", " ").split()[0]
                self.assertLessEqual(len(token), 6)
                self.assertEqual(float(token), threshold)
                self.assertEqual((changed/"tas").read_text(),
                                 (plain/"tas").read_text().replace("CHMAX=1.e-6", "CHMAX="+token))
            unsupported = root/"unsupported"
            with self.assertRaisesRegex(ValueError, "six-character field"):
                atmosphere_inputs(unsupported, prepared,
                                  dict(spec, temperature_convergence=1.23456e-7),
                                  table, abundance, masses, 4600, 5.65)
            self.assertFalse(unsupported.exists())

    def test_energy_balance_depth_is_explicit_and_inside_the_atmosphere(self):
        from generate_nongrey_grid import energy_balance_tau_division
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root/"data").mkdir()
            table = root/"table.bin"
            table.write_bytes(b"synthetic table")
            spec = dict(alpha=1.9, atmosphere_frequencies=5000, depths=100,
                        tau_top=1e-7, tau=100., tau_bottom=1000.)
            prepared = dict(data_sha256="synthetic", synple=str(root))
            abundance, masses = composition(.999, 0, [0]*5)
            plain, changed = root/"plain", root/"changed"
            atmosphere_inputs(plain, prepared, spec, table, abundance, masses, 4800, 6.)
            atmosphere_inputs(changed, prepared, dict(spec, energy_balance_tau_division=.0001),
                              table, abundance, masses, 4800, 6.)
            self.assertEqual((changed/"tas").read_text(),
                             (plain/"tas").read_text().replace("TAUDIV=0.01", "TAUDIV=0.0001"))
            for name in ["fort.5", "ember-masses.dat", "opacity.sha256"]:
                self.assertEqual((plain/name).read_bytes(), (changed/name).read_bytes())
            self.assertNotEqual(input_fingerprint("same-binary", plain),
                                input_fingerprint("same-binary", changed))
            for value in [0., spec["tau_top"], spec["tau"], math.inf, math.nan]:
                with self.subTest(value=value), self.assertRaisesRegex(ValueError, "transition depth"):
                    energy_balance_tau_division(dict(spec, energy_balance_tau_division=value))

    def test_temperature_iteration_stays_below_independent_correction_limit(self):
        from generate_nongrey_grid import temperature_convergence
        self.assertEqual(temperature_convergence({}), 1e-6)
        self.assertEqual(temperature_convergence(dict(temperature_convergence=1e-8)), 1e-8)
        self.assertEqual(temperature_convergence(dict(temperature_convergence=5e-6)), 5e-6)
        for value in [0., -1e-8, 1e-5, math.inf, math.nan]:
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "temperature convergence"):
                temperature_convergence(dict(temperature_convergence=value))
        with self.assertRaisesRegex(ValueError, "six-character field"):
            temperature_convergence(dict(temperature_convergence=1.23456e-7))

    def test_single_composition_uses_four_independent_opacity_workers(self):
        barrier = threading.Barrier(4, timeout=5)
        lock = threading.Lock()
        active = 0
        peak = 0
        seen = []

        def source(job):
            nonlocal active, peak
            with lock:
                active += 1
                peak = max(peak, active)
                seen.append(job)
            barrier.wait()
            with lock:
                active -= 1
            return job

        for planes in [["one"], ["one", "two"]]:
            seen.clear()
            expected = [(p, i, 1000 + 100*i) for p in planes
                        for i in reversed(range(8))]
            result = run_opacity_jobs(planes, list(range(1000, 1800, 100)), source, 8)
            self.assertEqual(result, expected)
            self.assertCountEqual(seen, expected)
            self.assertEqual(peak, 4)
            self.assertEqual(active, 0)

    def test_failed_isotherm_cannot_be_returned_as_a_complete_opacity_family(self):
        def source(job):
            if job[1] == 1:
                raise RuntimeError("source row failed validation")
            return job

        with self.assertRaisesRegex(RuntimeError, "failed validation"):
            run_opacity_jobs(["one"], [1000, 2000, 3000], source, 4)

    def test_trial_column_scaling_retains_temperature_and_hydrostatic_relation(self):
        n=24
        mass=[10**(-5+6*i/(n-1)) for i in range(n)]
        rows=[[3000+100*i,1e10*(i+1),1e-7*(i+1),1e13*(i+1)] for i in range(n)]
        seed=f'{n} -4\n'+'\n'.join(map(str,mass))+'\n'+'\n'.join(' '.join(map(str,r)) for r in rows)+'\n'
        self.assertEqual(scale_initial_column(seed,1),seed)
        trial=continuation_structure(seed,3000,5.4,3200,5.7)
        scaled=scale_initial_column(trial,2).split();original=trial.split()
        for i in range(n):
            self.assertEqual(float(scaled[2+i]),2*float(original[2+i]))
            for j in range(4):
                k=2+n+4*i+j
                self.assertEqual(float(scaled[k]),float(original[k])*(1 if j==0 else 2))
            # n*T and g*m increase by the same factor at the target gravity.
            k=2+n+4*i
            pressure_ratio=float(scaled[k])*float(scaled[k+3])/(float(original[k])*float(original[k+3]))
            self.assertAlmostEqual(pressure_ratio,float(scaled[2+i])/float(original[2+i]))
        for factor in [0,-1,math.inf,math.nan]:
            with self.assertRaisesRegex(ValueError,'finite and positive'):
                scale_initial_column(seed,factor)

    def test_plan_cancellation_checks_work_and_original_plan(self):
        from run_nongrey_plan import check_cancellation, PlanCancelled
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            check_cancellation(work, 'original-plan')
            record = {'work': str(work.resolve()), 'plan_sha256': 'original-plan',
                      'reason': 'Stop future source launches; preserve running models.'}
            path = work/'cancellation.json'
            path.write_text(json.dumps(record))
            with self.assertRaisesRegex(PlanCancelled, 'Stop future'):
                check_cancellation(work, 'original-plan')
            with self.assertRaisesRegex(ValueError, 'does not identify'):
                check_cancellation(work, 'different-plan')
            record['work'] = str(work.resolve()/'another-plan')
            path.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, 'does not identify'):
                check_cancellation(work, 'original-plan')

    def test_cancelled_plan_never_launches_a_source_process(self):
        script = Path(__file__).resolve().parents[1]/'scripts/run_nongrey_plan.py'
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            plan = work/'input.json'
            plan.write_text(json.dumps({'requests': [{'name': 'unstarted-model'}]}))
            (work/'cancellation.json').write_text(json.dumps({
                'work': str(work.resolve()),
                'plan_sha256': hashlib.sha256(plan.read_bytes()).hexdigest(),
                'reason': 'Preserve the completed source models.'}))
            result = subprocess.run([sys.executable, '-B', str(script), str(plan), str(work)],
                                    capture_output=True, text=True, timeout=20)
            self.assertNotEqual(result.returncode, 0)
            manifest = json.loads((work/'manifest.json').read_text())
            self.assertEqual(manifest['models'][0]['status'], 'cancelled')
            self.assertEqual((work/'plan.json').read_bytes(), plan.read_bytes())
            self.assertFalse((work/'unstarted-model.log').exists())

    def test_seed_bottom_uses_measured_optical_depth(self):
        n=40
        tau=[10**(-5+9*i/(n-1)) for i in range(n)]
        mass=[2*t for t in tau]
        rows=[[3000+t**.25, 1e12, 1e-5, 1e15] for t in tau]
        seed=f'{n} -4\n'+'\n'.join(map(str,mass))+'\n'+'\n'.join(' '.join(map(str,r)) for r in rows)+'\n'
        log='FINAL MODEL ATMOSPHERE\n'+'\n'.join(f'{i+1} {m} {t} 3000 1e12 1e-5 1e7 -3 .5 .5 1' for i,(m,t) in enumerate(zip(mass,tau)))
        result=truncate_initial_structure(seed,log,1000)
        words=result.split(); count=int(words[0])
        self.assertLess(count,n)
        self.assertEqual(float(words[2]),mass[0])
        self.assertAlmostEqual(float(words[1+count]),2000,places=9)
        self.assertEqual(int(resample_initial_structure(result,n).split()[0]),n)
        with self.assertRaisesRegex(ValueError,'strictly inside'):
            truncate_initial_structure(seed,log,1e5)
        with self.assertRaisesRegex(ValueError,'mass grid'):
            truncate_initial_structure(seed,log.replace(f'1 {mass[0]} ', '1 99 '),1000)

        # The native fort.7 checkpoint rounds mass to E15.6; run.log keeps
        # the full source value. This is serialization, not a new mass grid.
        printed=[format(m,'.6e') for m in mass]
        rounded=f'{n} -4\n'+'\n'.join(printed)+'\n'+'\n'.join(' '.join(map(str,r)) for r in rows)+'\n'
        truncated=truncate_initial_structure(rounded,log,1000).split()
        self.assertEqual(int(truncated[0]),count)
        self.assertAlmostEqual(float(truncated[1+count]),2000,delta=.001)
        damaged=rounded.replace(printed[1],format(float(printed[1])*1.000002,'.6e'),1)
        with self.assertRaisesRegex(ValueError,'mass grid'):
            truncate_initial_structure(damaged,log,1000)

    def test_seed_top_preserves_deeper_mesh(self):
        n=60; tau=[10**(-9+12*i/(n-1)) for i in range(n)]
        mass=[2*t for t in tau]; rows=[[3000+t**.25,1e12,1e-5,1e15]for t in tau]
        seed=f'{n} -4\n'+'\n'.join(map(str,mass))+'\n'+'\n'.join(' '.join(map(str,r))for r in rows)+'\n'
        log='FINAL MODEL ATMOSPHERE\n'+'\n'.join(f'{i+1} {m} {t} 3000 1e12 1e-5 1e7 -3 .5 .5 1'for i,(m,t)in enumerate(zip(mass,tau)))
        words=truncate_initial_structure(seed,log,top_tau=1e-6).split()
        cut=next(i for i,t in enumerate(tau)if t>=1e-6); count=n-cut
        self.assertEqual(int(words[0]),count)
        self.assertEqual(list(map(float,words[2:2+count])),mass[cut:])
        self.assertEqual(list(map(float,words[2+count:])),[x for row in rows[cut:]for x in row])
        both=truncate_initial_structure(seed,log,100.,top_tau=1e-6).split(); count=int(both[0])
        self.assertEqual(float(both[2]),mass[cut]); self.assertAlmostEqual(float(both[1+count]),200.)
        for top in (0.,-1.,float('nan'),1000.,100.):
            with self.assertRaises(ValueError):truncate_initial_structure(seed,log,100.,top_tau=top)
        with self.assertRaisesRegex(ValueError,'mass grid'):
            truncate_initial_structure(seed,log.replace(f'1 {mass[0]} ','1 99 '),top_tau=1e-6)
        with self.assertRaisesRegex(ValueError,'twenty'):
            truncate_initial_structure(seed,log,top_tau=50.)

    def test_warm_composition_extension_preserves_complete_old_cells(self):
        old_axes = [[.3, .7], [0, .12], [2600, 2800, 3000, 3200], [4.9, 5.15, 5.4]]
        old = set(itertools.product(*old_axes))
        added = set(itertools.product([.2], old_axes[1], [3000, 3200], old_axes[3]))
        axes = [[.2, .3, .7], *old_axes[1:]]
        cells = complete_cells(axes, old | added)
        self.assertEqual([cell for cell in cells if cell[0][0] == .3], complete_cells(old_axes, old))
        extension = [cell for cell in cells if cell[0][0] == .2]
        self.assertEqual(len(extension), 2)
        self.assertTrue(all(cell[2] == [3000, 3200] for cell in extension))
        # A lone accepted model cannot stand in for an interpolation cell.
        lone = complete_cells(axes, old | {(.2, 0, 3200, 5.15)})
        self.assertEqual(len(lone), len(complete_cells(old_axes, old)))

    def setUp(self):
        self.opacity={"temperature_K":[1000,10000],"density_g_cm3":[1e-13,1e-2]}
        self.flux=5.670374419e-5*2800**4/(4*math.pi)
        # Deliberately synthetic hydrostatic/flux records, not a source grid.
        self.log=f"TOTAL SURFACE FLUX {self.flux:.17e}\nFINAL MODEL ATMOSPHERE\n"
        self.correction=""
        for i in range(1,31):
            tau=10**(-4+7*(i-1)/29)
            self.log+=f"{i} {tau:.17e} {tau:.17e} 4000 1e12 1e-6 {tau*1e5:.17e} -3 .8 .2 1\n"
            self.correction+=f"8 {i} 1e-7 0 0 0 1e-7 0 0\n"
            self.log+=f"EMBER CHEMICAL DENSITY: {i} 4000 {tau*1e5:.17e} 1e-6 1e-6\n"

    def parse(self,log=None,correction=None,opacity=None):
        return source_state(log if log is not None else self.log,
            correction if correction is not None else self.correction,2800,5,
            opacity if opacity is not None else self.opacity)

    def test_matching_depth(self):
        s=self.parse()
        self.assertAlmostEqual(s["T"],4000,places=8)
        self.assertAlmostEqual(s["Pgas"],1e7,places=5)
        self.assertLess(s["tau_bracket"][0],100)
        self.assertGreater(s["tau_bracket"][1],100)

    def test_legacy_depth_metadata_preserves_exact_recorded_diagnostics(self):
        current=self.parse()
        added={'optical_depth_range','column_mass_range'}
        legacy={k:v for k,v in current.items() if k not in added}
        self.assertTrue(source_diagnostics_match(current,current))
        self.assertTrue(source_diagnostics_match(legacy,current))
        for saved in [current,legacy]:
            for name in ['T','Pgas','flux_error','depths','temperature_range']:
                missing=dict(saved);missing.pop(name)
                self.assertFalse(source_diagnostics_match(missing,current))
            altered=dict(saved);altered['T']*=1.0000000001
            self.assertFalse(source_diagnostics_match(altered,current))
            self.assertFalse(source_diagnostics_match({**saved,'unknown':1},current))
        partial=dict(current);partial.pop('optical_depth_range')
        self.assertFalse(source_diagnostics_match(partial,current))
        changed=dict(current);changed['optical_depth_range']=[1e-6,100]
        self.assertFalse(source_diagnostics_match(changed,current))

    def test_rejects_wrong_flux_even_with_small_correction(self):
        with self.assertRaisesRegex(ValueError,"flux error"):
            self.parse(log=self.log.replace(".8 .2 1\n",".82 .2 1.02\n"))

    def test_rejects_inconsistent_molecular_density(self):
        with self.assertRaisesRegex(ValueError,"density closure"):
            self.parse(log=self.log.replace("1e-6 1e-6\n","1e-6 1.5e-6\n"))
        with self.assertRaisesRegex(ValueError,"missing independent chemical"):
            self.parse(log="\n".join(line for line in self.log.splitlines() if not line.startswith("EMBER CHEMICAL")))

    def test_rejects_failed_fortran_stop(self):
        for message in ['NOT CONVERGE IN RUSSEL','STOP partf; temp<1000 K',
                        'STOP EMBER DEPLETION FAILED',' **** STOP in SOLVE after ITER 12','h2minus:Stop']:
            # Even an earlier apparent final profile cannot make a later
            # fatal source message acceptable.
            with self.subTest(message=message),self.assertRaisesRegex(ValueError,"failed"):
                self.parse(log=self.log+message+'\n')

    def test_rejects_missing_convergence(self):
        with self.assertRaises(ValueError): self.parse(correction="")
        with self.assertRaises(ValueError): self.parse(correction=self.correction[:-32])

    def test_rejects_source_support(self):
        with self.assertRaisesRegex(ValueError,"source support"):
            self.parse(opacity={"temperature_K":[5000,10000],"density_g_cm3":[1e-13,1e-2]})

    def test_rejects_wrong_teff_and_gravity(self):
        with self.assertRaisesRegex(ValueError,"Teff"):
            self.parse(log=self.log.replace(f"{self.flux:.17e}",f"{self.flux*1.1:.17e}"))
        with self.assertRaisesRegex(ValueError,"gravity"):
            source_state(self.log,self.correction,2800,5.1,self.opacity)

    def test_baryon_and_helium_number_normalization(self):
        for x,y in [(.7,0),(.55,.1),(.45,.12)]:
            a,w=composition(x,y,[.003,.0001,.001,.01,.0059])
            mu=1.66053906660e-24; hm=1.67333e-24
            self.assertAlmostEqual(sum(v*m for v,m in zip(a,w))*hm*x/mu,1.,places=13)
            self.assertAlmostEqual(a[1]*x,y/3+(1-x-y-.02)/4,places=14)
            self.assertAlmostEqual(a[0]*x,x,places=14)
            self.assertAlmostEqual(a[1]*w[1]*hm*x/mu,1-x-.02,places=14)

    def test_initial_structure_preserves_hydrostatic_power_law(self):
        mass = [10**(-4+5*i/19) for i in range(20)]
        seed = "20 -4\n" + "\n".join(str(v) for v in mass) + "\n"
        seed += "\n".join(f"{3000*m**.2} {1e10*m} {1e-5*m**.8} {1e15*m**.8}" for m in mass)
        refined = list(map(float,resample_initial_structure(seed,61).split()))
        self.assertEqual(refined[:2],[61,-4])
        self.assertAlmostEqual(refined[2],mass[0])
        self.assertAlmostEqual(refined[62],mass[-1])
        for i,m in enumerate(refined[2:63]):
            t,ne,rho,nt = refined[63+4*i:67+4*i]
            self.assertAlmostEqual(t/(3000*m**.2),1,places=13)
            self.assertAlmostEqual(rho*t/(.03*m),1,places=13)
        with self.assertRaises(ValueError): resample_initial_structure(seed,10)
        with self.assertRaises(ValueError): resample_initial_structure(seed.replace("20 -4","20 -3"),61)
        fine=resample_initial_structure(seed,61)
        with self.assertRaises(ValueError):resample_initial_structure(fine,20)
        coarsened=list(map(float,resample_initial_structure(fine,20,allow_coarsen=True).split()))
        self.assertEqual(coarsened[:2],[20,-4])
        for i,m in enumerate(coarsened[2:22]):
            t,ne,rho,nt=coarsened[22+4*i:26+4*i]
            self.assertAlmostEqual(m/mass[i],1,places=13)
            self.assertAlmostEqual(t/(3000*m**.2),1,places=13)
            self.assertAlmostEqual(rho*t/(.03*m),1,places=13)
        with self.assertRaises(ValueError):resample_initial_structure(fine,10,allow_coarsen=True)

    def test_continuation_preserves_hydrostatic_pressure(self):
        mass=[10**(-4+5*i/19) for i in range(20)]
        boltzmann=1.380649e-16
        rows=[]
        for m in mass:
            temperature=3000*m**.2
            particles=10**5.15*m/(boltzmann*temperature)
            rows.append([temperature,.001*particles,2e-24*particles,particles])
        seed='20 -4\n'+'\n'.join(map(str,mass))+'\n'+'\n'.join(' '.join(map(str,row)) for row in rows)
        for teff,logg in [(2600,5.4),(3200,4.9)]:
            result=list(map(float,continuation_structure(seed,2800,5.15,teff,logg).split()))
            for i,m in enumerate(result[2:22]):
                t,ne,rho,nt=result[22+4*i:26+4*i]
                self.assertAlmostEqual(nt*boltzmann*t/(10**logg*m),1,places=13)
                self.assertAlmostEqual(nt*t/(rows[i][3]*rows[i][0]),1,places=13)
                self.assertAlmostEqual(rho/nt,2e-24,delta=1e-37)
                self.assertAlmostEqual(ne/nt,.001,places=14)

    def test_source_cache_survives_tlusty_empty_output_unit(self):
        # TLUSTY creates fort.2 even though this full-table recipe does not
        # read it. A finished expensive atmosphere must remain reusable.
        with tempfile.TemporaryDirectory() as temporary:
            d = Path(temporary)
            (d/"fort.5").write_text("2800 5.15\nT T\n")
            (d/"tas").write_text("IFRYB=1\n")
            before = input_fingerprint("executable-digest",d)
            (d/"fort.2").touch()
            self.assertEqual(input_fingerprint("executable-digest",d),before)
            (d/"tas").write_text("IFRYB=0\n")
            self.assertNotEqual(input_fingerprint("executable-digest",d),before)

    def test_input_abundances_and_mass_loading(self):
        spec={"metals":[.003,.0001,.001,.01,.0059],"alpha":1.9,
              "atmosphere_frequencies":20000,"depths":200,"tau_top":1e-7,"tau_bottom":1000}
        a,w=composition(.55,.1,spec["metals"])
        inputs={"atmosphere_input":"2800 5\nT T\n 'tas'\n0\n92\n"+
                "".join(f"{1 if i < 30 else 0} {v:.17e} 0\n" for i,v in enumerate(a[:92]))+"0 0 0 -1 0 0 ' ' ' '\n",
                "element_masses":"\n".join(str(v) for v in w),
                "parameters":"IOPTAB=-1,IFRYB=1,IFMOL=1,TMOLIM=10000,IDLST=0,IFRAYL=1\n"
                             "HMIX0=1.9,IFRSET=20000,ND=200,TAUFIR=1e-7,TAULAS=1000,TAUDIV=.01,CHMAX=1e-6,ILGDER=1,NITER=200,DPSILT=1.03,DERT=.001"}
        marker=f"EMBER ELEMENT MASSES: {w[0]} {w[1]}\nEMBER MOLECULAR EQUILIBRIUM TOLERANCE: 1e-8\n"
        source_inputs(inputs,spec,.55,.1,2800,5,marker)
        precise = dict(spec, molecular_equilibrium_tolerance=1e-11)
        precise_marker = marker.replace("TOLERANCE: 1e-8", "TOLERANCE: 1e-11")
        source_inputs(inputs,precise,.55,.1,2800,5,precise_marker)
        for declared, actual in [(spec, precise_marker), (precise, marker),
                                 (precise, precise_marker+marker)]:
            with self.assertRaisesRegex(ValueError,"molecular equilibrium tolerance mismatch"):
                source_inputs(inputs,declared,.55,.1,2800,5,actual)
        for value in [0, -1e-11, 1e-7, math.nan, math.inf]:
            with self.assertRaisesRegex(ValueError,"invalid source molecular equilibrium tolerance"):
                source_inputs(inputs,dict(spec,molecular_equilibrium_tolerance=value),
                              .55,.1,2800,5,marker)
        tightened = {**inputs, "parameters":inputs["parameters"].replace("CHMAX=1e-6", "CHMAX=1e-8")}
        source_inputs(tightened,dict(spec,temperature_convergence=1e-8),.55,.1,2800,5,marker)
        bounded = {**inputs, "parameters":inputs["parameters"].replace("CHMAX=1e-6", "CHMAX=5e-6")}
        source_inputs(bounded,dict(spec,temperature_convergence=5e-6),.55,.1,2800,5,marker)
        with self.assertRaisesRegex(ValueError,"settings mismatch"):
            source_inputs(inputs,dict(spec,temperature_convergence=1e-8),.55,.1,2800,5,marker)
        balanced = {**inputs, "parameters":inputs["parameters"].replace("TAUDIV=.01", "TAUDIV=.0001")}
        with self.assertRaisesRegex(ValueError,"settings mismatch"):
            source_inputs(balanced,spec,.55,.1,2800,5,marker)
        source_inputs(balanced,dict(spec,energy_balance_tau_division=.0001),.55,.1,2800,5,marker)
        with self.assertRaisesRegex(ValueError,"settings mismatch"):
            source_inputs(inputs,dict(spec,energy_balance_tau_division=.0001),.55,.1,2800,5,marker)
        with self.assertRaisesRegex(ValueError,"settings mismatch"):
            source_inputs(inputs,dict(spec,newton_relaxation=.3),.55,.1,2800,5,marker)
        source_inputs({**inputs,"parameters":inputs["parameters"]+",ORELAX=.3"},
                      dict(spec,newton_relaxation=.3),.55,.1,2800,5,marker)
        with self.assertRaisesRegex(ValueError,"settings mismatch"):
            source_inputs({**inputs,"parameters":inputs["parameters"]+",ORELAX=1"},
                          dict(spec,newton_relaxation=.3),.55,.1,2800,5,marker)
        with self.assertRaisesRegex(ValueError,"abundance mismatch"):
            source_inputs(inputs,spec,.55,0,2800,5,marker)
        with self.assertRaisesRegex(ValueError,"did not load"):
            source_inputs(inputs,spec,.55,.1,2800,5,"")
        with self.assertRaisesRegex(ValueError,"molecular equilibrium"):
            source_inputs(inputs,spec,.55,.1,2800,5,marker.replace("TOLERANCE: 1e-8","TOLERANCE: .001"))
        with self.assertRaisesRegex(ValueError,"settings mismatch"):
            source_inputs({**inputs,"parameters":inputs["parameters"].replace("HMIX0=1.9","HMIX0=1.0")},
                          spec,.55,.1,2800,5,marker)

    def test_binary_opacity_records_and_composition(self):
        a,_=composition(.55,.1,[.003,.0001,.001,.01,.0059])
        def rec(value):
            n=struct.pack("<i",len(value));return n+value+n
        blob=b"".join(rec(struct.pack("<4sdd",b"X   ",v,v)) for v in a[:92])
        blob+=rec(struct.pack("<id",1,10000))+rec(struct.pack("<10i",*([1]*10)))
        blob+=rec(struct.pack("<3i",2,2,2))
        for values in [[math.log(2000),math.log(6000)],[-20.,-10.],[20.,21.,22.,23.]]:
            blob+=rec(struct.pack(f"<{len(values)}d",*values))
        for f in [1e15,1e14]:
            blob+=rec(struct.pack("<d",f))
            for r in range(2):blob+=rec(struct.pack("<2f",r+.5,r+1.5))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"opacity.bin";p.write_bytes(blob)
            table=validate_table(p,a,[2000,6000],[math.exp(-20),math.exp(-10)])
            self.assertEqual(table["shape"],(2,2,2))
            self.assertEqual(list(table["log_opacity"]),[.5,1.5,1.5,2.5]*2)
            with self.assertRaisesRegex(ValueError,"composition mismatch"):
                validate_table(p,[v*2 for v in a],[2000,6000],[math.exp(-20),math.exp(-10)])
            from assemble_nongrey_grid import validate_shared_table
            from unittest.mock import patch
            checksum=hashlib.sha256(blob).hexdigest()
            arguments=[p,checksum,a,[2000,6000],[math.exp(-20),math.exp(-10)]]
            with patch('assemble_nongrey_grid.validate_table',wraps=validate_table) as reader:
                validate_shared_table(*arguments)
                alias=Path(d)/'same-opacity';alias.symlink_to(p)
                validate_shared_table(alias,*arguments[1:])
                self.assertEqual(reader.call_count,1)
                with self.assertRaisesRegex(ValueError,'composition mismatch'):
                    validate_shared_table(p,checksum,[v*2 for v in a],*arguments[3:])
                with self.assertRaisesRegex(ValueError,'support mismatch'):
                    validate_shared_table(p,checksum,a,[2000,7000],arguments[-1])
            p.write_bytes(blob[:-1])
            with self.assertRaisesRegex(ValueError,'checksum mismatch'):
                validate_shared_table(*arguments)
            with self.assertRaisesRegex(ValueError,'truncated'):
                validate_shared_table(p,hashlib.sha256(p.read_bytes()).hexdigest(),*arguments[2:])
            with self.assertRaisesRegex(ValueError,"truncated"):read_table(p)
            p.write_bytes(blob+b"x")
            with self.assertRaisesRegex(ValueError,"trailing"):read_table(p)

    def test_isotherms_preserve_source_cell_order(self):
        def rec(fmt,*values):
            data=struct.pack("<"+fmt,*values);m=struct.pack("<i",len(data));return m+data+m
        for count in [3,33]:
            with self.subTest(temperature_rows=count),tempfile.TemporaryDirectory() as d:
                paths=[]
                for i in range(count):
                    t=2000*4**(i/(count-1))
                    blob=b"".join(rec("4sdd",b"X   ",1.,1.) for _ in range(92))
                    blob+=rec("id",1,10000)+rec("10i",*([1]*10))+rec("3i",2,1,2)
                    blob+=rec("d",math.log(t))+rec("2d",-20,-10)+rec("2d",20+i,25+i)
                    for k,f in enumerate([1e15,1e14]):
                        blob+=rec("d",f)
                        for j in range(2):blob+=rec("f",100*k+10*j+i)
                    p=Path(d)/str(i);p.write_bytes(blob);paths.append(p)
                output=Path(d)/"combined"
                merge_isotherms(paths,output)
                table=read_table(output)
                self.assertEqual(table["shape"],(2,2,count))
                self.assertEqual(list(table["log_opacity"]),[base+i for base in [0,10,100,110] for i in range(count)])
                self.assertEqual(table["log_electron_density"],tuple(v+i for i in range(count) for v in [20,25]))
                with self.assertRaisesRegex(ValueError,"unordered"):
                    merge_isotherms(list(reversed(paths)),output)

    def test_archived_reference_flux_convergence(self):
        root=Path(__file__).resolve().parents[1]/"data/atmosphere/sources/nongrey_validation"
        opacity={"temperature_K":[3000,15000],"density_g_cm3":[1e-12,1e-5]}
        for depths in [70,100,200]:
            log=gzip.decompress((root/f"benchmark-{depths}-run.log.gz").read_bytes()).decode()
            changes=gzip.decompress((root/f"benchmark-{depths}-fort.9.gz").read_bytes()).decode()
            if depths < 200:
                with self.assertRaisesRegex(ValueError,"flux error"):
                    source_state(log,changes,5500,4.5,opacity,require_chemical_closure=False)
            else:
                s=source_state(log,changes,5500,4.5,opacity,require_chemical_closure=False)
                self.assertAlmostEqual(s["T"],9889.132721137456,places=7)
                self.assertAlmostEqual(s["Pgas"],263084.78873140743,places=6)

    def test_archived_electron_grid_matches_independent_state_log(self):
        root=Path(__file__).resolve().parents[1]/"data/atmosphere/sources/nongrey_validation"
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"opacity"
            for mode in ["sampled","dense"]:
                p.write_bytes(gzip.decompress((root/f"opacity-{mode}-fort.63.gz").read_bytes()))
                table=read_table(p);nt=table["shape"][2];nr=table["shape"][1]
                states=gzip.decompress((root/f"opacity-{mode}-fort.29.gz").read_bytes()).decode()
                count=0
                for line in states.splitlines():
                    r=line.split()
                    if len(r)!=10 or not r[0].isdigit():continue
                    it,ir=int(r[0])-1,int(r[1])-1
                    actual=math.exp(table["log_electron_density"][it*nr+ir])
                    self.assertLess(abs(actual/float(r[4])-1),.006)
                    count+=1
                self.assertEqual(count,nt*nr)

    def test_archived_cool_molecular_density_closure(self):
        root=Path(__file__).resolve().parents[1]/"data/atmosphere/sources/nongrey_validation"
        log=gzip.decompress((root/"chemical-pilot-100-run.log.gz").read_bytes()).decode()
        changes=gzip.decompress((root/"chemical-pilot-100-fort.9.gz").read_bytes()).decode()
        s=source_state(log,changes,2800,5.15,
                       {"temperature_K":[1400,15000],"density_g_cm3":[1e-13,1e-3]})
        self.assertLess(s["chemical_density_error"],2e-8)
        self.assertAlmostEqual(s["T"],4031.7664955028854,places=7)
        self.assertAlmostEqual(s["source_density"],9.691865923733929e-5,places=12)

    def test_archived_density_rollback_control_is_rejected(self):
        root=Path(__file__).resolve().parents[1]/"data/atmosphere/sources/nongrey_validation"
        log=gzip.decompress((root/"rollback-control-run.log.gz").read_bytes()).decode()
        changes=gzip.decompress((root/"rollback-control-fort.9.gz").read_bytes()).decode()
        with self.assertRaisesRegex(ValueError,"molecular density closure"):
            source_state(log,changes,2800,5.15,
                         {"temperature_K":[1400,15000],"density_g_cm3":[1e-13,1e-3]})

    def test_installed_physical_family_reimports_without_missing_cells(self):
        data=Path(__file__).resolve().parents[1]/"data/atmosphere"
        root=data/"sources/nongrey_gs98_z020"
        manifest=json.loads((root/"manifest.json").read_text())
        for name,expected in manifest["archive_files_sha256"].items():
            self.assertEqual(hashlib.sha256((root/name).read_bytes()).hexdigest(),expected,name)
        for name,expected in manifest["provenance"]["patches"].items():
            self.assert_patch_digest(data/"sources"/name, expected)
        with tempfile.TemporaryDirectory() as temporary:
            d=Path(temporary)
            records=import_grid(root/"manifest.json",d/"grid.dat")
            self.assertEqual(len(records),48)
            actual=shlex.split((d/"grid.dat").read_text())
            expected=shlex.split((data/"nongrey_gs98_z020_tau100.dat").read_text())
            self.assertEqual(len(actual),len(expected))
            for a,b in zip(actual,expected):
                try:
                    x,y=float(a),float(b)
                except ValueError:
                    self.assertEqual(a,b)
                else:
                    # Allow last-bit libm differences across build hosts.
                    self.assertTrue(math.isclose(x,y,rel_tol=2e-14,abs_tol=1e-15),(a,b))
            for record in manifest["models"]:
                for kind in ["log","convergence","atmosphere_input","element_masses","parameters","initial_structure"]:
                    if kind in record: record[kind]=str(root/record[kind])
            manifest["models"].pop()
            (d/"manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,"incomplete physical source grid"):
                import_grid(d/"manifest.json",d/"bad.dat")
            manifest["models"][0]["log_sha256"]="0"*64
            (d/"manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,"output checksum mismatch"):
                import_grid(d/"manifest.json",d/"bad.dat")

    def test_archived_completed_structure_can_seed_a_new_source_run(self):
        from generate_nongrey_grid import completed_initial_structure
        root=Path(__file__).resolve().parents[1]/"data/atmosphere/sources/nongrey_gs98_z020"
        spec=json.loads((root/"specification.json").read_text())
        model=root/"plane-003/model-003-002"
        with tempfile.TemporaryDirectory() as temporary:
            d=Path(temporary)
            for name in ["completed.json","fort.7","fort.9","run.log"]:
                (d/(name+".gz")).write_bytes((model/(name+".gz")).read_bytes())
            initial=completed_initial_structure(d,spec,3200,5.4,300)
            self.assertEqual(initial,gzip.decompress((model/"fort.7.gz").read_bytes()).decode())
            # A saved structure with a valid gzip container but altered bytes
            # cannot borrow the accepted atmosphere's receipt.
            (d/"fort.7.gz").write_bytes(gzip.compress((initial+"\n").encode()))
            with self.assertRaisesRegex(ValueError,"receipt checksum mismatch"):
                completed_initial_structure(d,spec,3200,5.4,300)

    def test_archived_native_convection_jacobian(self):
        root=Path(__file__).resolve().parents[1]
        report=json.loads((root/"docs/results/nongrey_convection_jacobian.json").read_text())
        sources=root/"data/atmosphere/sources"
        for name,expected in report["files"].items():
            self.assertEqual(hashlib.sha256((sources/"nongrey_validation"/name).read_bytes()).hexdigest(),expected)
        archived_patch=gzip.decompress((sources/"nongrey_validation"/report["corrected_patch_archive"]).read_bytes())
        self.assertEqual(hashlib.sha256(archived_patch).hexdigest(),
                         report["corrected_patch_sha256"])
        self.assertEqual(hashlib.sha256((root/"scripts/nongrey_convection_probe.f").read_bytes()).hexdigest(),
                         report["corrected"]["driver_sha256"])
        for mode in ["original","corrected"]:
            text=gzip.decompress((sources/"nongrey_validation"/f"convection-jacobian-{mode}-output.txt.gz").read_bytes()).decode()
            rows=[list(map(float,line.split())) for line in text.splitlines()]
            self.assertEqual(rows,report[mode]["rows"])
            self.assertEqual(len(rows),9)
            for case,column,coefficient,finite_difference,error in rows:
                if mode=="corrected":
                    self.assertLess(abs(coefficient/finite_difference-1),2e-6)
                elif column==2:
                    self.assertLess(coefficient*finite_difference,0)

    def test_archived_native_opacity_derivative_and_state_restoration(self):
        root=Path(__file__).resolve().parents[1]
        report=json.loads((root/"docs/results/nongrey_opacity_derivative.json").read_text())
        sources=root/"data/atmosphere/sources/nongrey_validation"
        self.assert_patch_digest(sources.parent/"tlusty208-ember.patch",
                                 report["corrected_patch_sha256"])
        for name,expected in report["files"].items():
            self.assertEqual(hashlib.sha256((sources/name).read_bytes()).hexdigest(),expected)
        self.assertEqual(hashlib.sha256((root/"scripts/nongrey_opacity_derivative_probe.f").read_bytes()).hexdigest(),
                         report["corrected"]["driver_sha256"])
        for mode in ["original","corrected"]:
            output=gzip.decompress((sources/f"opacity-derivative-{mode}-output.txt.gz").read_bytes()).decode()
            rows=[list(map(float,line.split())) for line in output.splitlines()]
            self.assertEqual(rows,report[mode]["rows"])
            self.assertEqual(len(rows),6)
            for frequency,depth,derivative,temperature,density,mass,pressure in rows:
                self.assertTrue(all(math.isfinite(v) for v in [derivative,temperature,density,mass,pressure]))
                if mode=="corrected":
                    self.assertLess(max(abs(v) for v in [derivative,temperature,density,mass,pressure]),1e-10)
                else:
                    self.assertAlmostEqual(derivative,(depth+1)/4-1,places=12)
                    self.assertGreater(abs(mass),.001)
                    self.assertGreater(abs(pressure),.001)

    def test_interrupted_checkpoint_is_only_a_verified_initial_guess(self):
        from generate_nongrey_grid import checkpoint_initial_structure
        root=Path(__file__).resolve().parents[1]
        report=json.loads((root/"docs/results/nongrey_initial_checkpoint.json").read_text())
        sources=root/"data/atmosphere/sources/nongrey_validation"
        with tempfile.TemporaryDirectory() as temporary:
            d=Path(temporary)
            for name,expected in report["files"].items():
                raw=(sources/name).read_bytes()
                self.assertEqual(hashlib.sha256(raw).hexdigest(),expected)
                (d/name.removeprefix("initial-checkpoint-").removesuffix(".gz")).write_bytes(gzip.decompress(raw))
            text=checkpoint_initial_structure(d,.7,.12,3200,5.4)
            self.assertEqual(int(text.split()[0]),100)
            with self.assertRaisesRegex(ValueError,"no final atmosphere"):
                source_state((d/"run.log").read_text(),(d/"fort.9").read_text(),3200,5.4,
                             {"temperature_K":[1000,15000],"density_g_cm3":[1e-13,1e-3]})
            with self.assertRaisesRegex(ValueError,"grid label"):
                checkpoint_initial_structure(d,.45,.12,3200,5.4)
            (d/"fort.7").write_text(text+"\n")
            with self.assertRaisesRegex(ValueError,"checksum"):
                checkpoint_initial_structure(d,.7,.12,3200,5.4)


if __name__ == "__main__": unittest.main()
