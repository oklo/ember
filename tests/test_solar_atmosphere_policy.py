from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts'
sys.path[:0]=[str(SCRIPTS),str(SCRIPTS/'solar_grid')]
from boundary_replay import change_family
from ember_convection_seed import split_regions
from release import thermal_turns
from warm_boundary import thickness
from ember_campaign_recovery import tick,save


class AtmospherePolicyTest(unittest.TestCase):
    def test_saved_split_signatures(self):
        self.assertEqual(split_regions(150,list(range(164,186)),[149],209),[[149,162],[185,207]])
        self.assertEqual(split_regions(180,list(range(181,190)),list(range(153,158)),209),[[153,157],[189,207]])
        self.assertEqual(split_regions(151,[],[148,150],205),[[148,148],[150,203]])
        self.assertIsNone(split_regions(150,[155,158],[149],209))

    def test_remove_interior_knot_without_touching_endpoints(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            profile=d/'layer'
            profile.write_text('header\nscope\ntau 3\ncomposition\nlog_teff 3 2 3 4\nlog_g 2 4 5\ndata\na\nb\nc\nd\ne\nf\n')
            family=d/'family';family.write_text(f'file "{profile}"\n')
            change_family(family,d/'trial',1000)
            rows=(d/'trial/layer').read_text().splitlines()
            self.assertEqual(rows[4],'log_teff 2 2 4')
            self.assertEqual(rows[7:],['a','b','e','f'])

    def test_turns_retain_temperature_and_age(self):
        points=[dict(years=i,Teff_K=t) for i,t in enumerate([10,8,6,6.02,7,9,8])]
        self.assertEqual(thermal_turns(points),[dict(kind='minimum',years=2,Teff_K=6),dict(kind='maximum',years=5,Teff_K=9)])

    def test_skip_requires_live_prerequisites_and_releases_family(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            base=dict(trigger='working',command=['false'],outputs=[],wall_seconds=1)
            jobs=[dict(base,name='prepare'),dict(base,name='mid',after=['prepare'],
                  interpolate_when=['decision.json','interpolate_midpoint']),dict(base,name='family',after=['mid'])]
            save(root/'recovery-plan.json',jobs)
            save(root/'decision.json',{'interpolate_midpoint':True})
            def states():return tick(root,{6:('working',None)},0,False)[1]
            self.assertEqual(states()['mid'],'waiting')
            parent=root/'recovery/prepare';parent.mkdir(parents=True)
            save(parent/'result.json',{'status':'passed'})
            self.assertEqual(states()['mid'],'interpolated')
            self.assertEqual(states()['family'],'ready')

    def test_failed_optional_calibration_allows_measured_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            base=dict(trigger='working',command=['false'],outputs=[],wall_seconds=1)
            jobs=[dict(base,name='replay'),dict(base,name='decide',after_attempt=['replay'])]
            save(root/'recovery-plan.json',jobs)
            parent=root/'recovery/replay';parent.mkdir(parents=True)
            save(parent/'result.json',{'status':'failed'})
            self.assertEqual(tick(root,{6:('working',None)},0,False)[1]['decide'],'ready')

    def test_thickness_at_constant_density(self):
        data='FINAL MODEL ATMOSPHERE\n'+'\n'.join(f'{i} {m} {m} 1000 1 2 1 1 1 1 1' for i,m in enumerate([.1,1,2,4],1))
        self.assertAlmostEqual(thickness(data),(3-2/3)/2)
