"""Low-Z imports must retain distinct native planes and original source cells."""
from contextlib import redirect_stdout
import hashlib
import io
from pathlib import Path
import shlex
import sys
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile, ZIP_DEFLATED

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import import_aesopus_mixtures as importer


class LowMetalImport(unittest.TestCase):
    def test_zero_and_submilliz_sources_are_distinct_and_unmodified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); archive=root/'native.zip'; output=root/'output'
            metals=[0., .0001, .001]
            xs=[0., .1, .2, .35, .5, .7, .8, .9, .95]
            temperatures=list(range(2000,3701,10))+list(range(3720,4501,20))
            densities=list(range(-8000,6001,200))
            expected={z:f'{-1-i/10:.6f}' for i,z in enumerate(metals)}
            with ZipFile(archive,'w',ZIP_DEFLATED) as source:
                for z in metals:
                    for x in xs:
                        for region in ['lowT','highT','highR']:
                            rr=[r for r in densities if (r>1000)==(region=='highR')]
                            tt=[t for t in temperatures if region=='highR' or (t<=3700)==(region=='lowT')]
                            header=(f'# Calculations performed with AESOPUS 2.1\n'
                                    f'# TABLE X= {x:.6f} Zref= {z:.6f} Z= {z:.6f}\n'
                                    f'# log10(R) range: nre= {len(rr)} values from {rr[0]/1000:g} to {rr[-1]/1000:g}\n'
                                    f'# in steps of 0.2\n')
                            rows=''.join(f'{t/1000:.3f} '+ ' '.join([expected[z]]*len(rr))+'\n' for t in tt)
                            name=f'03-GS98/GS98_a0.0_OPALZ_{region}/aesopus2.0_gasbroad_GS98_Z{z:.6f}_X{x:g}.tab'
                            source.writestr(name,header+rows)
            digest=hashlib.sha256(archive.read_bytes()).hexdigest()
            argv=['import',str(archive),str(output),'--metals',*map(str,metals)]
            with patch.object(importer,'SOURCE_SHA256',digest), patch.object(sys,'argv',argv), redirect_stdout(io.StringIO()):
                importer.main()
            lines=(output/'aesopus21_gs98_mixture.dat').read_text().splitlines()
            planes=[shlex.split(s) for s in lines[1:]]
            self.assertEqual(len({name for _,name in planes}),3)
            self.assertEqual(planes[0][1],'aesopus21_gs98_mixture_z000.dat')
            self.assertEqual(planes[2][1],'aesopus21_gs98_mixture_z001.dat')
            for z_string,name in planes:
                z=float(z_string); content=(output/name).read_text().splitlines()
                self.assertEqual(content[0].split()[:3],['9','211','71'])
                for plane in range(len(xs)):
                    offset=3+plane*(len(temperatures)+1)
                    self.assertEqual(float(content[offset].split()[1]),z)
                    cells=[v for row in content[offset+1:offset+1+len(temperatures)] for v in row.split()]
                    self.assertEqual(len(cells),211*71)
                    self.assertEqual(set(cells),{expected[z]})


if __name__=='__main__':
    unittest.main()
