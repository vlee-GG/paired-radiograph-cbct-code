"""Reproduce five-view independent retrieval over the released fixed candidate sets.
NCC helpers derived from correspondence_check.py (MIT); see LICENSE.
"""
from pathlib import Path
import argparse,csv,json,concurrent.futures,collections,time,shutil
import numpy as np
from scipy import fft as scipy_fft
import five_view_core as ref

def score_mask(model,angle):
    vertical={60:62,75:50,90:27,105:50,120:62}[int(angle)]
    return (np.abs(ref.f0)[:,None]<=vertical)&(np.abs(ref.f1)[None,:]<=20)

def score_model(model,query_indices,ids,candidates,workers,outpath):
    if outpath.exists():return np.load(outpath,allow_pickle=False)
    grads=np.load(CACHE/'new_view_gradients.npy',mmap_mode='r');dims=np.load(CACHE/'ct_reference_dims.npy')
    reference=np.load(CACHE/'ct_reference_pixels.npy',mmap_mode='r')
    needed=collections.defaultdict(list)
    for local,qi in enumerate(query_indices):
        for slot,ci in enumerate(candidates[qi]):needed[int(ci)].append((local,qi,slot))
    scores=np.full((len(query_indices),4,20),np.nan,float)
    checkpoint=outpath.with_suffix('.checkpoint.npz')
    finished=set()
    if checkpoint.exists():
        with np.load(checkpoint,allow_pickle=False) as z:
            assert np.array_equal(z['query_indices'],query_indices)
            assert str(z['model_name'])==model['name']
            scores[:]=z['scores'];finished=set(map(int,z['completed_candidates']))
    uses_reference=model['axis'] is None or model.get('projection')=='reference'
    projection_indices=None
    def one(ci):
        h,w=dims[ci];uses=needed[ci];all_results=[]
        for ai in range(4):
            pix=reference[ci,:h,:w]
            vf=ref.variant_ffts(pix);values=np.empty(len(uses));mask=score_mask(model,ref.NEWANGLES[ai])
            for start in range(0,len(uses),16):
                part=uses[start:start+16];g=np.stack([grads[qi,ai] for local,qi,slot in part])
                a=g-g.mean(axis=(-2,-1),keepdims=True);nr=(a*a).sum(axis=(-2,-1)).astype(float)
                if np.any(nr<=0):raise ValueError('Zero-gradient image')
                fr=np.fft.rfft2(a,s=ref.SHAPE,axes=(-2,-1));best=np.full(len(part),-np.inf)
                for fd,nd in vf:
                    if nd<=0:raise ValueError('Zero-gradient projection')
                    cc=scipy_fft.irfft2(fr*np.conj(fd),s=ref.SHAPE,axes=(-2,-1),workers=1)
                    best=np.maximum(best,cc[:,mask].max(axis=1)/np.sqrt(nr*nd))
                values[start:start+len(part)]=best
            all_results.append(values)
        return ci,np.array(all_results)
    t=time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for ci,values in pool.map(one,[x for x in sorted(needed) if x not in finished]):
            for col,(local,qi,slot) in enumerate(needed[ci]):scores[local,:,slot]=values[:,col]
            finished.add(ci);count=len(finished)
            if count%100==0 or count==len(needed):
                temp=checkpoint.with_name(checkpoint.stem+'.tmp.npz')
                np.savez(temp,scores=scores,query_indices=np.array(query_indices),model_name=np.array(model['name']),completed_candidates=np.array(sorted(finished)))
                temp.replace(checkpoint)
                print('Model',model['name'],'candidates',count,'/',len(needed),'seconds',round(time.monotonic()-t),flush=True)
    assert np.isfinite(scores).all();np.save(outpath,scores);return scores
def main():
    global CACHE
    ap=argparse.ArgumentParser();ap.add_argument('--dataset',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=2);ap.add_argument('--cache',type=Path);ap.add_argument('--queries',type=Path,help='Optional sample_id CSV for a reproduction spot check');args=ap.parse_args()
    ref.ROOT=args.dataset.resolve();ref.OUT=args.out.resolve()
    if ref.OUT.is_relative_to(ref.ROOT):raise ValueError('Output must be outside dataset')
    ref.OUT.mkdir(parents=True,exist_ok=True);CACHE=args.cache.resolve() if args.cache else ref.OUT/'cache';ref.CACHE=CACHE
    material=ref.ROOT/'validation/multiview_retrieval'
    # New metadata/documentation files do not alter the frozen candidate bank.
    shutil.copyfile(material/'candidate_sets.csv',ref.OUT/'candidate_sets.csv')
    fingerprints={n:ref.sha(ref.ROOT/n) for n in ('metadata/metadata.csv','metadata/files.csv','metadata/partition_units.csv','metadata/splits/example_5fold_by_partition_unit.csv','SHA256SUMS.txt')}
    ref.savejson(ref.OUT/'frozen_protocol.json',dict(dataset_fingerprints=fingerprints))
    protocol,meta,ids,index,candidates,files=ref.load_inputs()
    baseline=ref.reuse_90(ids,candidates,index)
    if not args.cache:ref.prepare_features(meta,ids,files,args.workers)
    query=list(range(len(ids))) if not args.queries else [index[x['sample_id']] for x in ref.readcsv(args.queries)]
    model=dict(name='angle_shift_axis1_factor1',axis=1,factor=1.,projection='reference')
    new=score_model(model,query,ids,candidates,args.workers,ref.OUT/'scores_oblique.npy')
    five=np.stack([new[:,0],new[:,1],baseline[query],new[:,2],new[:,3]],axis=1)
    np.savez_compressed(ref.OUT/'candidate_scores.npz',ids=np.array(ids)[query],angles=np.array(ref.ANGLES),candidate_ids=np.array(ids)[candidates[query]],scores=five)
    ranks=1+(five[:,:,1:]>five[:,:,0,None]).sum(axis=2)
    with (ref.OUT/'summary.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['nominal_angle_deg','records','top1_percent','top3_percent','mean_rank'])
        for ai,angle in enumerate(ref.ANGLES):w.writerow([angle,len(query),100*(ranks[:,ai]==1).mean(),100*(ranks[:,ai]<=3).mean(),ranks[:,ai].mean()])
    print('Completed',len(query),'queries across all five angles',flush=True)

if __name__=='__main__':main()
