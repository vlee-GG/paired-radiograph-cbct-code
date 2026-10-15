"""Numerical helpers for five-view independent retrieval (MIT; see LICENSE)."""
from pathlib import Path
import os
import argparse,collections,concurrent.futures,csv,hashlib,json,time,os,gc
import numpy as np
import nibabel as nib
import pydicom
from PIL import Image
from scipy import fft as scipy_fft

ROOT=Path(os.environ.get('PAIRED_DATASET','physionet-paired-radiograph-cbct')).resolve()
OUT=Path(os.environ.get('PAIRED_OUT','multiview_results')).resolve()
if OUT.is_relative_to(ROOT):
    raise ValueError('Output must be outside the input dataset')
CACHE=OUT/'cache'
ANGLES=(60,75,90,105,120);NEWANGLES=(60,75,105,120)
RSIZE=(114,139);SHAPE=(278,228)
VARIANTS=((0,1.),(-6,1.),(6,1.),(-3,1.),(3,1.),(0,1.08),(-6,1.08),(6,1.08))
f0=np.fft.fftfreq(SHAPE[0],1./SHAPE[0]).astype(int)
f1=np.fft.fftfreq(SHAPE[1],1./SHAPE[1]).astype(int)
MASK=(np.abs(f0)[:,None]<=27)&(np.abs(f1)[None,:]<=20)

def readcsv(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def savejson(p,obj):
    tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(obj,indent=2),encoding='utf-8');tmp.replace(p)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def win(a,lo=1.,hi=99.5):
    l,h=np.percentile(a,lo),np.percentile(a,hi)
    return np.clip((a-l)/max(h-l,1e-6),0,1)

def grad_mag(x):
    gx=np.zeros_like(x);gy=np.zeros_like(x)
    gx[:,1:-1]=x[:,2:]-x[:,:-2];gy[1:-1,:]=x[2:,:]-x[:-2,:]
    g=np.sqrt(gx*gx+gy*gy);return g/(g.max()+1e-9)

def fft_e(g):
    a=g-g.mean();return np.fft.rfft2(a,s=SHAPE),float((a*a).sum())

def radio_gradient(path):
    im=pydicom.dcmread(path).pixel_array.astype(np.float32)
    im=Image.fromarray((win(im)*255).astype(np.uint8)).resize(RSIZE,Image.Resampling.BILINEAR)
    return grad_mag(np.asarray(im,np.float32)/255.)

def variant_ffts(pixels):
    img0=Image.fromarray(pixels);out=[]
    for rot,mag in VARIANTS:
        im=img0
        if mag>1:
            w,h=im.size;cw,ch=int(round(w/mag)),int(round(h/mag))
            im=im.crop(((w-cw)//2,(h-ch)//2,(w+cw)//2,(h+ch)//2))
        if rot:im=im.rotate(rot,resample=Image.Resampling.BILINEAR,fillcolor=0)
        g=grad_mag(np.asarray(im.resize(RSIZE,Image.Resampling.BILINEAR),np.float32)/255.)
        out.append(fft_e(g))
    return out

def load_inputs():
    protocol=json.loads((OUT/'frozen_protocol.json').read_text())
    for n,h in protocol['dataset_fingerprints'].items():
        if sha(ROOT/n)!=h:raise ValueError('Dataset metadata changed after input freezing: '+n)
    meta={r['sample_id']:r for r in readcsv(ROOT/'metadata/metadata.csv')};ids=sorted(meta);index={s:i for i,s in enumerate(ids)}
    candidates=np.full((len(ids),20),-1,np.int32)
    for r in readcsv(OUT/'candidate_sets.csv'):candidates[index[r['query_id']],int(r['candidate_index'])]=index[r['candidate_id']]
    assert np.all(candidates>=0) and np.array_equal(candidates[:,0],np.arange(len(ids)))
    files={r['relative_path']:r for r in readcsv(ROOT/'metadata/files.csv')}
    return protocol,meta,ids,index,candidates,files

def base(meta,sid):
    m=meta[sid];return ROOT/'data'/m['species']/m['data_group_id']/m['source_volume_id']/sid

def reuse_90(ids,candidates,index):
    old=Path(os.environ.get('PAIRED_BASELINE90',str(Path(__file__).resolve().parent/'reference'/'baseline_90.npz')))
    z=np.load(old,allow_pickle=False);assert list(z['ids'])==ids
    score=np.empty((len(ids),20),float);score[:,0]=z['matched']
    neg=z['mismatched'].reshape(len(ids),19)
    for qi in range(len(ids)):
        order=np.argsort(candidates[qi,1:])
        score[qi,1:][order]=neg[qi]
    published={r['sample_id']:r for r in readcsv(ROOT/'validation/correspondence_check.csv')}
    ranks=1+(score[:,1:]>score[:,0,None]).sum(axis=1)
    assert np.array_equal(ranks,z['ranks'])
    for qi,sid in enumerate(ids):
        r=published[sid];top=1+int(np.argmax(score[qi,1:]));bestid=ids[candidates[qi,top]]
        assert int(r['rank'])==ranks[qi]
        assert abs(float(r['score_matched'])-score[qi,0])<=.000051
        assert abs(float(r['score_best_other'])-score[qi,top])<=.000051
        assert r['best_other_id']==bestid,(sid,bestid,r['best_other_id'])
    np.save(OUT/'scores_a090_reused.npy',score)
    savejson(OUT/'baseline_90_reuse.json',dict(status='verified_and_reused',records=len(ids),candidate_scores=score.size,source_npz_sha256=sha(old),published_csv_sha256=sha(ROOT/'validation/correspondence_check.csv'),mapping='Historical distractor score order is sorted candidate identifier, as in the public candidate-first loop; mapped to the frozen sampled candidate slots.',checks=['all identifiers','all ranks','all matched and best-other rounded scores','all best-other identifiers'],full_90_scoring_repeated=False))
    print('Reused all 90-degree results; all published ranks/scores/best-other IDs agree',flush=True)
    return score

def prepare_features(meta,ids,files,workers):
    CACHE.mkdir(parents=True,exist_ok=True)
    complete=CACHE/'features_complete.json'
    if complete.exists():
        sig=json.loads(complete.read_text());assert sig['input_protocol_sha256']==sha(OUT/'frozen_protocol.json')
        print('Reusing completed local descriptor cache',flush=True);return
    grads=np.lib.format.open_memmap(CACHE/'new_view_gradients.npy',mode='w+',dtype=np.float32,shape=(len(ids),4,139,114))
    projections=np.lib.format.open_memmap(CACHE/'ct_reference_pixels.npy',mode='w+',dtype=np.uint8,shape=(len(ids),120,100))
    dims=np.zeros((len(ids),2),np.int32)
    def one(q):
        sid=ids[q];b=base(meta,sid);p=b/(sid+'__cbct_crop.nii.gz')
        assert sha(p)==files[p.relative_to(ROOT).as_posix()]['sha256']
        a=np.asanyarray(nib.load(p).dataobj).astype(np.float32)
        pix=(win(a.mean(axis=2)[10:110,20:140].T)*255).astype(np.uint8)
        gg=np.stack([radio_gradient(b/(sid+f'__xray__a{ang:03d}.dcm')) for ang in NEWANGLES])
        return q,pix,gg
    t=time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for count,(q,pix,gg) in enumerate(pool.map(one,range(len(ids))),1):
            h,w=pix.shape;assert h<=120 and w<=100
            projections[q,:h,:w]=pix;dims[q]=h,w;grads[q]=gg
            if count%100==0:print('Prepared CT reference + four views',count,'/',len(ids),'seconds',round(time.monotonic()-t),flush=True)
    grads.flush();projections.flush();np.save(CACHE/'ct_reference_dims.npy',dims)
    savejson(complete,dict(input_protocol_sha256=sha(OUT/'frozen_protocol.json'),records=len(ids),new_angles=list(NEWANGLES),elapsed_seconds=round(time.monotonic()-t,2)))
    del grads,projections;gc.collect()


"""Reproduce five-view independent retrieval over the released fixed candidate sets.
NCC helpers derived from correspondence_check.py (MIT); see LICENSE.
"""
from pathlib import Path
import argparse,csv,json,concurrent.futures,collections,time,shutil
import numpy as np
from scipy import fft as scipy_fft
import sys
ref=sys.modules[__name__]

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
