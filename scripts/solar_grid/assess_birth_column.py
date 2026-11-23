"""Assess whichever native calculation actually finished the shared column."""
from pathlib import Path
from atmosphere_repair import score_column
from coordinate import read,root,write
selection=read(root/'recovery/birth-column-warm-restart/selected-column.json')
directory=Path(selection['directory'])
result=score_column(directory)
write(root/'birth-atmosphere-2850/reviewed-column.json',dict(result,directory=str(directory)))
raise SystemExit(0 if result['passes'] else 1)
