#!/usr/bin/env python3
import copy
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from summarize_remnant_endpoints import summarize


class Endpoints(unittest.TestCase):
    def track(self):
        return {'converged':True,'columns':['age_yr','central_X','Teff_K','L_Lsun','nuclear_deposited_Lsun'],
                'history':[[0,.7,3000,1,.9],[10,1e-3,1500,1,.01],
                           [20,1e-6,400,1,1e-4],[30,1e-6,1200,1,.5],[40,1e-6,90,1,1e-5]]}

    def test_keeps_shell_reignition_and_temperature_recrossings(self):
        result=summarize(self.track())
        nuclear=result['nuclear_deposited_over_photon_luminosity'][0]['crossings']
        self.assertEqual([c['direction'] for c in nuclear],['down','up','down'])
        self.assertEqual(nuclear[0]['age_bracket_yr'],[0,10])
        temperature=result['temperature'][0]['crossings']
        self.assertEqual([c['direction'] for c in temperature],['down','up','down'])
        self.assertEqual(result['temperature'][-1]['crossings'][0]['age_bracket_yr'],[30,40])
        self.assertEqual(len(result['core_hydrogen'][0]['crossings']),1)

    def test_missing_legacy_nuclear_data_are_unknown(self):
        track=self.track();track['columns']=track['columns'][:-1]
        track['history']=[r[:-1] for r in track['history']]
        self.assertIsNone(summarize(track)['nuclear_deposited_over_photon_luminosity'])

    def test_continuation_already_below_is_not_an_exhaustion_age(self):
        track=self.track();track['history']=track['history'][2:]
        result=summarize(track)['core_hydrogen'][0]
        self.assertTrue(result['already_below_at_segment_start'])
        self.assertEqual(result['crossings'],[])

    def test_unordered_or_nonphysical_history_rejected(self):
        for value,index in [(-1,0),(float('nan'),1),(0,2),(-1,3)]:
            track=copy.deepcopy(self.track());track['history'][1][index]=value
            with self.assertRaises(ValueError):summarize(track)


if __name__=='__main__':unittest.main()
