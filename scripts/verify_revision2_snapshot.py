"""Verify archived bytes and denominators, not the validity of every inference."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

repo=Path(__file__).resolve().parents[1]
manifest=json.loads((repo/'provenance/revision2_sources.json').read_text())
for item in manifest['files']:
    path=repo/item['destination']
    if hashlib.sha256(path.read_bytes()).hexdigest()!=item['destination_sha256']:
        raise SystemExit(f'Checksum mismatch: {item["destination"]}')
def read(relative):
    with (repo/'results/revision2'/relative).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))
supp='supplementary_tables_S0-S5/'
s0=read(supp+'S0_per_run.csv')
assert len(s0)==2601
filtered=[r for r in s0 if r['feature_method']!='Variance']
assert Counter(r['result_set'] for r in filtered)=={'ncurve':1852,'abl':117,'hart':52,'baronfix':94}
s0b=read(supp+'S0b_rev2_per_run.csv')
assert len(s0b)==2097
assert Counter((r['result_set'],r['status']) for r in s0b)=={('abl2','ok'):234,('reps','ok'):1553,('krange','ok'):298,('krange','skipped_k_below_2'):12}
assert len(read(supp+'S5_missing_configurations.csv'))==20
assert Counter(r['result_set'] for r in read(supp+'S5b_rev2_missing.csv'))=={'reps':7,'krange':2}
for path,n in [('07_simulation/E8_assumption.csv',300),('07_simulation/sim_complete_per_panel.csv',1263),('08_markers/marker_coverage_per_dataset.csv',3120),('08_markers/marker_coverage_tests.csv',768)]:
    assert len(read(path))==n,path
for name in ['S1_ARI_per_run.csv','S2_FM_per_run.csv','S3_Jaccard_per_run.csv','S4_Purity_per_run.csv']:
    assert len([r for r in read(supp+name) if r['feature_method']!='Variance'])==1852,name
assert hashlib.md5((repo/'inst/extdata/dip_null_tables.rds').read_bytes()).hexdigest()=='29406eeaabe9806321450b034074c868'
print(f'PASS: {len(manifest["files"])} source/result hashes and all documented denominators.')
print('For numerical inference checks, also run scripts/verify_revision2_statistics.py.')
