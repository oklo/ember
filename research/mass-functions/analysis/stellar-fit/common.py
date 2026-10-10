from pathlib import Path
import json,sys
import numpy as np
W=Path(__file__).resolve().parent;V2=W.parent/'cluster-dissolution-oct9-v2';V1=W.parent/'cluster-dissolution-oct9-v1';C=W.parent/'cluster-mass-functions-oct9-v1'
sys.path.insert(0,str(V2))
import analysis_functions as old
CAL={(r['cluster'],r['band']):dict(coefficients=r['coefficients'],xmin=r['logL_min'],xmax=r['logL_max'],scatter=max(.12,r['robust_scatter_mag'])) for r in json.loads((V2/'photometric_calibration.json').read_text())}
ANGLES={'hyades':np.array([0,10,20,40,80,180.000001]),'pleiades':np.array([0,1,2,3,5,10,30,90,180.000001])}
MASS=np.array([.03,.08,.10,.125,.16,.20])
STAR_BINS=[2,3,4];BD_BIN=0
KMS=1.022712165045695
