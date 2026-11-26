import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/solar_grid'))
from adjudicate import record
from release import decision, orders, routing, save
from shepherd import current


class ReleaseTest(unittest.TestCase):
    def test_adoption_invalidated_by_changed_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            order = dict(id='surface',targets=['pilot'],restart_on_adoption=['next'])
            save(root/'work-orders.json',[order])
            save(root/'status.json',dict(recovery={'pilot':'passed'}))
            save(root/'assessment.json',{'passes':True})
            jobs = dict(pilot={'after':['column']},column={},next={})
            priorities,gates=routing(root,jobs)
            self.assertIn('column',priorities)
            self.assertIn('next',gates)
            record(root,'surface','adopted','reviewer','matches withheld column',['assessment.json'])
            self.assertEqual(decision(root,order),'adopted')
            self.assertNotIn('next',routing(root,jobs)[1])
            save(root/'assessment.json',{'passes':False})
            self.assertEqual(decision(root,order),'review stale')
            self.assertIn('next',routing(root,jobs)[1])

    def test_pending_target_cannot_be_adopted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            save(root/'work-orders.json',[dict(id='clock',targets=['crossing'])])
            save(root/'status.json',dict(recovery={'crossing':'running'}))
            with self.assertRaises(ValueError):
                record(root,'clock','adopted','reviewer','conditional',['status.json'])

    def test_stale_heartbeat_cannot_request_provider(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            save(root/'status.json',dict(hold=None,phase='working',updated_unix=0))
            self.assertFalse(current(root,dict(complete=0)))


if __name__ == '__main__':
    unittest.main()
