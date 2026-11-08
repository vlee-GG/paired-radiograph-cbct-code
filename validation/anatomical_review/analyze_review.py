"""Recompute the first-submission study from the released CSV tables."""
from pathlib import Path
import csv,json,argparse,collections
import numpy as np

def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig',newline='')))
def cluster_ci(rows,role,key):
    units=sorted({r[key] for r in rows});ix={v:i for i,v in enumerate(units)};arr=np.zeros((len(units),4))
    for r in rows:
        j=ix[r[key]];o=r['reference_class']=='original';arr[j,0 if o else 2]+=1
        arr[j,1 if o else 3]+=r[role+'_answer']==('yes' if o else 'no')
    rng=np.random.default_rng(20261020);out=[]
    for _ in range(40):
        x=arr[rng.integers(len(units),size=(500,len(units)))].sum(1);x=x[(x[:,0]>0)&(x[:,2]>0)]
        out.append(np.stack([(x[:,1]+x[:,3])/(x[:,0]+x[:,2]),x[:,1]/x[:,0],x[:,3]/x[:,2],.5*(x[:,1]/x[:,0]+x[:,3]/x[:,2])],1))
    z=np.concatenate(out)
    return {'clusters':len(units),'replicates':len(z),'seed':20261020,'intervals':{k:np.quantile(z[:,j],[.025,.975]).tolist() for j,k in enumerate(['overall','original_yes','substituted_no','balanced'])}}

def main(root):
    rows=read(root/'review_records.csv');p=read(root/'presentation_records.csv');saved=json.loads((root/'review_stats.json').read_text())
    assert len(rows)==300 and len({r['composite_id'] for r in rows})==300
    result={}
    for role in ['reader1','reader2']:
        s=saved['readers'][role]
        for cl in ['original','substituted']:
            cc=collections.Counter(r[role+'_answer'] for r in rows if r['reference_class']==cl)
            assert all(cc[a]==s['by_reference_class'][cl][a] for a in ['yes','no','uncertain'])
        correct=sum(r[role+'_answer']==('yes' if r['reference_class']=='original' else 'no') for r in rows)
        assert correct==s['label_agreement_n']
        for key,target in [('partition_unit_id','partition_cluster_CI'),('source_volume_id','source_cluster_CI')]:
            actual=cluster_ci(rows,role,key)
            for k,v in actual['intervals'].items():assert np.allclose(v,s[target]['intervals'][k]),(role,target,k)
        seen={};agree=0;repeats=0
        for r in p:
            if r['reader']!=role or r['block']!='formal':continue
            cid=r['composite_id']
            if cid in seen:repeats+=1;agree+=seen[cid]==r['answer']
            else:seen[cid]=r['answer']
        assert repeats==15 and agree==s['hidden_repeat_agreement_n']
        result[role]={'label_agreement_n':correct,'n':300,'repeat_agreement':agree,'repeat_pairs':repeats}
    cats=['yes','no','uncertain'];matrix=np.zeros((3,3),int)
    for r in rows:matrix[cats.index(r['reader1_answer']),cats.index(r['reader2_answer'])]+=1
    po=np.trace(matrix)/300;p1=matrix.sum(1)/300;p2=matrix.sum(0)/300;pe=(p1*p2).sum();pool=(p1+p2)/2;pc=(pool*(1-pool)).sum()/2
    assert matrix.tolist()==saved['between_readers']['matrix']
    assert np.isclose((po-pe)/(1-pe),saved['between_readers']['kappa'])
    assert np.isclose((po-pc)/(1-pc),saved['between_readers']['AC1_3category'])
    assert sum(r['readers_disagree']=='True' for r in rows)==89
    follow=read(root/'author_review_records.csv');fstats=saved['author_followup']
    scope={r['composite_id']:r for r in rows if any(r[k]!=('yes' if r['reference_class']=='original' else 'no') for k in ['reader1_answer','reader2_answer'])}
    assert {r['composite_id'] for r in rows if r['author_review_in_scope']=='True'}==set(scope)
    assert {cid for cid,r in scope.items() if r['readers_disagree']!='True'}=={'c0065','c0280'}
    assert len(follow)==91 and len({r['composite_id'] for r in follow})==91 and {r['composite_id'] for r in follow}==set(scope)
    basis=collections.Counter();outcomes={}
    for r in follow:
        fixed=scope[r['composite_id']]
        for k in ['reference_class','xray_record_id','ct_crop_record_id','reader1_answer','reader2_answer']:assert r[k]==fixed[k]
        assert r['status']=='complete' and r['reviewer'] and r['review_date'] and r['basis_codes']
        codes=r['basis_codes'].split(';');assert set(codes)<=set(fstats['basis_options']);basis.update(codes)
        assert r['conclusion'] in ['corresponding','different','insufficient']
        assert not(r['conclusion']=='different' and codes==['limited']),r['composite_id']
        if r['conclusion']=='corresponding':assert {'bone','joint'}&set(codes)
    for cl in ['original','substituted']:outcomes[cl]=dict(collections.Counter(r['conclusion'] for r in follow if r['reference_class']==cl))
    assert outcomes==fstats['outcomes_by_reference_class'] and dict(basis)==fstats['basis_counts']
    assert all(r['labels_visible']=='True' for r in follow)
    assert fstats['cases_reviewed']==91 and fstats['original_scores_replaced'] is False
    concordant=sum(r['conclusion']==('corresponding' if r['reference_class']=='original' else 'different') for r in follow)
    assert concordant==fstats['reference_label_agreement_n']==91
    assert fstats['reference_label_agreement_denominator']==len(follow)==91 and fstats['reference_label_agreement_percent']==100.0
    print(json.dumps({'verified':True,'readers':result,'between_reader_agreement':float(po),'disagreements':89,'author_followup_cases':len(follow),'shared_reference_deviations':2},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent);a=p.parse_args();main(a.root)
