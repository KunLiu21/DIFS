"""Create a separate run directory; never overwrite different existing files."""
import argparse
import hashlib
import shutil
from pathlib import Path

repo = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--workspace', type=Path, default=repo/'work/revision2')
args = parser.parse_args()
dest = args.workspace.resolve()
if dest == repo or repo.is_relative_to(dest):
    parser.error('Choose a separate workspace, not the repository or an ancestor')
copies = []
for group in ('benchmark','simulation','inst'):
    copies += [(p,dest/p.relative_to(repo)) for p in (repo/group).rglob('*') if p.is_file()]
for folder, name in [('01_main_grid','grid_corrected.csv'),('02_ablation','abl_corrected.csv')]:
    copies.append((repo/'results/revision2'/folder/name, dest/'results/diag'/name))
copies.append((repo/'inst/extdata/dip_null_tables.rds',dest/'results/dip_null_tables.rds'))
conflicts = [str(q) for p,q in copies if q.exists() and (not q.is_file() or p.read_bytes()!=q.read_bytes())]
if conflicts:
    parser.error('No files copied; different files already exist: '+', '.join(conflicts[:10]))
for p,q in copies:
    q.parent.mkdir(parents=True,exist_ok=True)
    if not q.exists(): shutil.copyfile(p,q)
for d in ('source','scenarios','results/abl2','results/reps','results/krange','results/markers','results/sim_complete'):
    (dest/d).mkdir(parents=True,exist_ok=True)
print(f'Workspace ready: {dest}; {len(copies)} source/support files matched or copied')
print('No datasets downloaded or experiments executed. Raw per-run RDS files are not bundled.')
