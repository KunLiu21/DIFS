"""Independently check corrected signed-rank summaries using subset-sum counts."""
import csv
import math
from collections import Counter
from pathlib import Path

root=Path(__file__).resolve().parents[1]/'results/revision2'
def read(path):
    with (root/path).open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
def number(x): return float('nan') if x in ('','NA') else float(x)
def check(values, marker=False):
    assert all(math.isfinite(x) for x in values)
    d=[round(x,12) for x in values]
    z=[x for x in d if x!=0]
    p=float('nan') if marker and len(z)<2 else 1.
    if z and not (marker and len(z)<2):
        order=sorted(range(len(z)),key=lambda i:abs(z[i])); ranks=[0]*len(z)
        i=0
        while i<len(z):
            j=i+1
            while j<len(z) and abs(z[order[j]])==abs(z[order[i]]): j+=1
            for k in order[i:j]: ranks[k]=i+1+j  # doubled midrank
            i=j
        dist=Counter({0:1})
        for r in ranks:
            nxt=Counter(dist)
            for k,v in dist.items(): nxt[k+r]+=v
            dist=nxt
        total=sum(ranks); obs=sum(r for r,x in zip(ranks,z) if x>0)
        p=sum(v for k,v in dist.items() if abs(2*k-total)>=abs(2*obs-total))/(2**len(z))
    return p,sum(x>0 for x in d),sum(x<0 for x in d)
assert check([1,2,3])[0]==.25
assert check([1e-16,1,2,3])==check([0,1,2,3])
assert check([1/11,-(2/11-1/11)])[0]==1
def eq(a,b): return (math.isnan(a) and math.isnan(b)) or abs(a-b)<1e-12

ind={'Kumar','Trapnell','Romanov','Lawlor','Koh','Darmanis','Muraro','Fletcher','Baron'}
arms={}
for row in read('02_ablation/abl_corrected.csv'):
    key=(row['data.name'],int(row['n']),row['fmethod']); assert key not in arms
    arms[key]=number(row['ARI'])
for row in read('supplementary_tables_S0-S5/S0b_rev2_per_run.csv'):
    if row['result_set']!='abl2' or row['status']!='ok': continue
    assert row['stage2only_fell_back']!='TRUE'
    key=(row['dataset'],int(row['n_requested']),row['arm']); assert key not in arms
    arms[key]=number(row['ARI'])
comparisons={
 'dip ranking over the gate':('DIFS_stage1only','GateOnly'),
 'HVG ranking over the gate':('HVG_stage1only','GateOnly'),
 'stage I: dip minus HVG':('DIFS_stage1only','HVG_stage1only'),
 '+ stage II (no MSE), dip':('DIFS_fixed','DIFS_stage1only'),
 '+ MSE ratio search, dip':('DIFS_fixratio','DIFS_fixed'),
 '+ MSE size sweep, dip':('DIFS','DIFS_fixratio'),
 '+ stage II (no MSE), HVG':('HVG_fixed','HVG_stage1only'),
 '+ MSE (size and ratio), HVG':('HVG_full','HVG_fixed'),
 'full pipeline: dip minus HVG':('DIFS','HVG_full'),
 'no-MSE pipeline: dip minus HVG':('DIFS_fixed','HVG_fixed'),
 'stage II only minus stage I only':('DIFS_stage2only','DIFS_stage1only'),
 'combined minus stage II only':('DIFS','DIFS_stage2only')}
ab=read('02_ablation/rev2_abl2_comparisons.csv')
assert len(ab)==72
for row in ab:
    a,b=comparisons[row['comparison']]; n=int(row['n'])
    ds=sorted(d for d,k,arm in arms if k==n and arm==a and (row['datasets']=='13' or d in ind))
    diffs=[arms[d,n,a]-arms[d,n,b] for d in ds]
    p,pos,neg=check(diffs)
    assert eq(p,number(row['p'])) and pos==int(row['ahead']) and len(ds)==int(row['n_ds']),row
    assert abs(sum(diffs)/len(diffs)-number(row['mean']))<1e-12,row

markers={}
for row in read('08_markers/marker_coverage_per_dataset.csv'):
    if row['budget'] in ('','NA'): continue
    key=(row['dataset'],row['list'],row['method'],int(row['budget']))
    assert key not in markers
    markers[key]=row
tests=read('08_markers/marker_coverage_tests.csv');assert len(tests)==768
for row in tests:
    L=row['list']; b=int(row['budget']);comp=row['vs'];metric=row['metric']
    ds=sorted(d for d,l,m,k in markers if l==L and m=='DIFS' and k==b and (row['datasets']=='all' or d in ind) and (d,L,comp,b) in markers)
    diffs=[number(markers[d,L,'DIFS',b][metric])-number(markers[d,L,comp,b][metric]) for d in ds]
    p,pos,neg=check(diffs,marker=True)
    assert eq(p,number(row['p_exact'])) and pos==int(row['wins']) and neg==int(row['losses']) and len(ds)==int(row['n']),row
    assert abs(sum(round(x,12) for x in diffs)/len(diffs)-number(row['mean_diff']))<1e-12,row
print('PASS: all 72 component and 768 marker comparisons reproduce from per-run/per-dataset observations.')
print('Checked means, paired sample counts, zero/tie-aware win/loss counts and exact p-values; not a clustering rerun.')
