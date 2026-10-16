"""Recompute all five view retrieval scores from current images on fixed candidates.

No stored 90-degree scores or pre-existing image-feature caches enter computation.
"""
from pathlib import Path
import argparse, collections, concurrent.futures, csv, gzip, hashlib, io, json, os, shutil, time
from datetime import datetime, timezone
import numpy as np
import nibabel as nib
import pydicom
from PIL import Image
from scipy import fft as scipy_fft
from scipy.stats import rankdata

ANGLES=(60,75,90,105,120)
SHAPE=(278,228);RSIZE=(114,139)
VARIANTS=((0,1.),(-6,1.),(6,1.),(-3,1.),(3,1.),(0,1.08),(-6,1.08),(6,1.08))
VERTICAL=(62,50,27,50,62)
f0=np.fft.fftfreq(SHAPE[0],1./SHAPE[0]).astype(int)
f1=np.fft.fftfreq(SHAPE[1],1./SHAPE[1]).astype(int)
MASKS=[(np.abs(f0)[:,None]<=v)&(np.abs(f1)[None,:]<=20) for v in VERTICAL]
ROOT=Path('../physionet-paired-radiograph-cbct')
OUT=Path('./results').resolve()
def sha(data):return hashlib.sha256(data).hexdigest()
def file_sha(p):return sha(Path(p).read_bytes())
def readcsv(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def savecsv(path,rr):
    with Path(path).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def savejson(path,data):
    p=Path(path);t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(data,indent=2,ensure_ascii=False),encoding='utf-8');t.replace(p)
def progress(stage,done,total,**kw):
    r=dict(stage=stage,completed=done,total=total,updated_at=datetime.now(timezone.utc).isoformat(),**kw)
    savejson(OUT/'progress.json',r);print(stage,done,'/',total,kw,flush=True)
def win(a):
    lo,hi=np.percentile(a,[1.,99.5]);return np.clip((a-lo)/max(hi-lo,1e-6),0,1)
def gradient(a):
    gx=np.zeros_like(a);gy=np.zeros_like(a)
    gx[:,1:-1]=a[:,2:]-a[:,:-2];gy[1:-1]=a[2:]-a[:-2]
    g=np.sqrt(gx*gx+gy*gy);return g/(g.max()+1e-9)
def radio_gradient(array):
    im=Image.fromarray((win(array)*255).astype(np.uint8)).resize(RSIZE,Image.Resampling.BILINEAR)
    return gradient(np.asarray(im,np.float32)/255.)
def variant_ffts(pix):
    image=Image.fromarray(pix);out=[]
    for rot,mag in VARIANTS:
        p=image
        if mag>1:
            w,h=p.size;cw,ch=int(round(w/mag)),int(round(h/mag));p=p.crop(((w-cw)//2,(h-ch)//2,(w+cw)//2,(h+ch)//2))
        if rot:p=p.rotate(rot,resample=Image.Resampling.BILINEAR,fillcolor=0)
        g=gradient(np.asarray(p.resize(RSIZE,Image.Resampling.BILINEAR),np.float32)/255.)
        a=g-g.mean();out.append((np.fft.rfft2(a,s=SHAPE),float((a*a).sum())))
    return out
def subset(m):return 'mouse' if m['species']=='mouse' else {'limb_unspecified':'limb','maxillary':'maxilla','mandibular':'mandible'}[m['anatomical_region']]
def base(m,sid):return ROOT/'data'/m['species']/m['data_group_id']/m['source_volume_id']/sid

def freeze():
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'cache').mkdir(exist_ok=True);(OUT/'input_snapshot').mkdir(exist_ok=True)
    names=['metadata/metadata.csv','metadata/files.csv','metadata/partition_units.csv','metadata/splits/example_5fold_by_partition_unit.csv','SHA256SUMS.txt','validation/multiview_retrieval/candidate_sets.csv']
    protocol_path=OUT/'protocol.json'
    fingerprints={n:file_sha(ROOT/n) for n in names}
    if protocol_path.exists():
        protocol=json.loads(protocol_path.read_text(encoding='utf-8'))
        if fingerprints!=protocol['dataset_fingerprints']:raise ValueError('Dataset changed since this run was frozen; use a new output directory.')
    else:
        for name in names:
            p=OUT/'input_snapshot'/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/name,p)
        protocol=dict(analysis='fresh_five_view_independent_full_release',created_at=datetime.now(timezone.utc).isoformat(),dataset=str(ROOT),dataset_fingerprints=fingerprints,angles=ANGLES,candidates_per_query=20,seed=20260920,working_shape=[139,114],fft_shape=SHAPE,ct_reference='float32 NIfTI crop mean along axis 2; [10:110,20:140].T; percentile 1-99.5',radiograph='decoded float32; percentile 1-99.5; uint8; PIL bilinear resize 114 columns by 139 rows; gradient magnitude',variants=VARIANTS,horizontal_shift=20,vertical_shift=VERTICAL,rank='1 + count(distractors strictly greater than recorded assignment)',parameter_selection='No tuning or new selection in this rerun; prior published configuration retained.',fresh_90_recomputation=True,prior_scores_reused_for_computation=False,old_image_features_reused=False,score_dtype='float64',scope='All release queries against the same fixed 20 candidates per record; not exhaustive all-candidate retrieval or calibrated registration.')
        savejson(protocol_path,protocol)
    meta={r['sample_id']:r for r in readcsv(OUT/'input_snapshot/metadata/metadata.csv')};ids=sorted(meta);index={s:i for i,s in enumerate(ids)}
    candidates=np.full((len(ids),20),-1,np.int32)
    for r in readcsv(OUT/'input_snapshot/validation/multiview_retrieval/candidate_sets.csv'):
        candidates[index[r['query_id']],int(r['candidate_index'])]=index[r['candidate_id']]
    assert len(ids)==2229 and np.all(candidates>=0) and np.array_equal(candidates[:,0],np.arange(len(ids)))
    for qi,c in enumerate(candidates):
        assert len(set(map(int,c)))==20
        assert all(subset(meta[ids[x]])==subset(meta[ids[qi]]) for x in c)
    files={r['relative_path']:r for r in readcsv(OUT/'input_snapshot/metadata/files.csv')}
    manifest={line.split(maxsplit=1)[1].lstrip('*').strip():line.split(maxsplit=1)[0] for line in (OUT/'input_snapshot/SHA256SUMS.txt').read_text(encoding='utf-8-sig').splitlines()}
    return protocol,meta,ids,candidates,files,manifest

def prepare_features(meta,ids,files,manifest,workers):
    cache=OUT/'cache';complete=cache/'features_complete.json'
    if complete.exists():
        assert json.loads(complete.read_text())['protocol_sha256']==file_sha(OUT/'protocol.json')
        print('Resume: features generated earlier in THIS frozen run.',flush=True);return
    grads=np.lib.format.open_memmap(cache/'radiograph_gradients.npy',mode='w+',dtype=np.float32,shape=(len(ids),5,139,114))
    ref=np.lib.format.open_memmap(cache/'ct_reference_pixels.npy',mode='w+',dtype=np.uint8,shape=(len(ids),120,100))
    dims=np.zeros((len(ids),2),np.int32);input_rows=[];qc_rows=[];t0=time.monotonic()
    def checked(path,role,sid,angle=''):
        before=path.stat();data=path.read_bytes();after=path.stat();h=sha(data);rel=path.relative_to(ROOT).as_posix()
        assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns),('File changed while reading',rel)
        expected=files[rel]
        if h!=expected['sha256'] or len(data)!=int(expected['bytes']) or h!=manifest[rel]:raise ValueError('Input hash/size differs from frozen metadata: '+rel)
        row=dict(sample_id=sid,angle=angle,role=role,relative_path=rel,bytes=len(data),sha256=h,mtime_ns=after.st_mtime_ns)
        return data,row
    def one(qi):
        sid=ids[qi];b=base(meta[sid],sid);entries=[];quality=[]
        raw,r=checked(b/(sid+'__cbct_crop.nii.gz'),'ct_crop',sid);entries.append(r)
        image=nib.Nifti1Image.from_bytes(gzip.decompress(raw));a=np.asanyarray(image.dataobj).astype(np.float32)
        pix=(win(a.mean(axis=2)[10:110,20:140].T)*255).astype(np.uint8)
        gg=[]
        for angle in ANGLES:
            raw,r=checked(b/(sid+f'__xray__a{angle:03d}.dcm'),'radiograph',sid,angle);entries.append(r)
            ds=pydicom.dcmread(io.BytesIO(raw));arr=ds.pixel_array.astype(np.float32)
            assert arr.shape==(834,688) and np.isfinite(arr).all() and np.ptp(arr)>0,(sid,angle)
            g=radio_gradient(arr);assert np.isfinite(g).all() and np.var(g)>0
            gg.append(g);quality.append(dict(sample_id=sid,angle=angle,rows=arr.shape[0],columns=arr.shape[1],decode_ok=True,finite=True,nonconstant=True,sha256_verified=True,bytes_verified=True,sop_uid=str(ds.SOPInstanceUID)))
        return qi,pix,np.stack(gg),entries,quality
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for count,(qi,pix,gg,entries,quality) in enumerate(pool.map(one,range(len(ids))),1):
            h,w=pix.shape;ref[qi,:h,:w]=pix;dims[qi]=h,w;grads[qi]=gg;input_rows.extend(entries);qc_rows.extend(quality)
            if count%100==0 or count==len(ids):progress('fresh_image_features',count,len(ids),seconds=round(time.monotonic()-t0,1))
    grads.flush();ref.flush();np.save(cache/'ct_reference_dims.npy',dims)
    savecsv(OUT/'input_file_hashes.csv',input_rows);savecsv(OUT/'radiograph_integrity.csv',qc_rows)
    savejson(complete,dict(protocol_sha256=file_sha(OUT/'protocol.json'),records=len(ids),image_files=len(input_rows),radiographs=len(qc_rows),elapsed_seconds=time.monotonic()-t0))

def score_all(ids,candidates,workers,fft_workers):
    cache=OUT/'cache';grads=np.load(cache/'radiograph_gradients.npy',mmap_mode='r');ref=np.load(cache/'ct_reference_pixels.npy',mmap_mode='r');dims=np.load(cache/'ct_reference_dims.npy')
    score_path=OUT/'scores_full.npy';done_path=OUT/'completed_candidates.npy'
    if score_path.exists() and done_path.exists():
        scores=np.load(score_path,mmap_mode='r+');done=np.load(done_path)
    else:
        scores=np.lib.format.open_memmap(score_path,mode='w+',dtype=np.float64,shape=(len(ids),5,20));scores[:]=np.nan;done=np.zeros(len(ids),bool)
    needed=[[] for _ in ids]
    for qi,row in enumerate(candidates):
        for slot,ci in enumerate(row):needed[int(ci)].append((qi,slot))
    def one(ci):
        h,w=map(int,dims[ci]);transforms=variant_ffts(ref[ci,:h,:w]);uses=needed[ci];values=np.empty((len(uses),5),np.float64)
        for start in range(0,len(uses),8):
            part=uses[start:start+8];query=[q for q,slot in part];g=np.array(grads[query],copy=True)
            a=g-g.mean(axis=(-2,-1),keepdims=True);norm=(a*a).sum(axis=(-2,-1)).astype(float)
            if np.any(norm<=0):raise ValueError('Zero-energy query')
            fr=np.fft.rfft2(a,s=SHAPE,axes=(-2,-1));best=np.full((len(part),5),-np.inf)
            for fd,energy in transforms:
                cc=scipy_fft.irfft2(fr*np.conj(fd),s=SHAPE,axes=(-2,-1),workers=fft_workers)
                for ai,mask in enumerate(MASKS):best[:,ai]=np.maximum(best[:,ai],cc[:,ai,mask].max(axis=1)/np.sqrt(norm[:,ai]*energy))
            values[start:start+len(part)]=best
        return ci,values
    t0=time.monotonic();initial=int(done.sum());count=initial
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(one,ci):ci for ci in range(len(ids)) if not done[ci]}
        for future in concurrent.futures.as_completed(futures):
            ci,values=future.result()
            for k,(qi,slot) in enumerate(needed[ci]):scores[qi,:,slot]=values[k]
            done[ci]=True;count+=1
            if count%25==0 or count==len(ids):
                scores.flush();np.save(done_path,done)
                progress('fresh_five_view_scores',count,len(ids),seconds=round(time.monotonic()-t0,1),resumed=initial)
    assert np.isfinite(scores).all() and np.all(done)
    return np.asarray(scores)

def verify_against_reference(ids,candidates,scores):
    # Independently repeat a small set through the prior angle-first numeric loop.
    cache=OUT/'cache';grads=np.load(cache/'radiograph_gradients.npy',mmap_mode='r');ref=np.load(cache/'ct_reference_pixels.npy',mmap_mode='r');dims=np.load(cache/'ct_reference_dims.npy')
    query=[ids.index(s) for s in ('mou_s001','mou_s033','mou_s062','mou_s111','rab_s0001','rab_s0351','rab_s0996','rab_s1001')]
    largest=0.;n=0
    for qi in query:
        for ci_slot,ci in enumerate(candidates[qi]):
            h,w=map(int,dims[ci]);vf=variant_ffts(ref[ci,:h,:w])
            for ai in range(5):
                g=np.array(grads[qi,ai]);a=g-g.mean();fr=np.fft.rfft2(a,s=SHAPE);nr=float((a*a).sum());best=-np.inf
                for fd,nd in vf:
                    cc=scipy_fft.irfft2(fr*np.conj(fd),s=SHAPE,workers=1)
                    best=max(best,float(cc[MASKS[ai]].max()/np.sqrt(nr*nd)))
                largest=max(largest,abs(best-float(scores[qi,ai,ci_slot])));n+=1
    assert largest<1e-12,largest
    savejson(OUT/'independent_loop_verification.json',dict(status='passed',queries=8,scores=n,max_absolute_difference=largest,method='Original angle-first single-score loop vs all-angle batched rerun; all scores freshly recomputed.'))

def finalize(meta,ids,candidates,scores):
    ranks=1+(scores[:,:,1:]>scores[:,:,0,None]).sum(axis=2)
    np.savez_compressed(OUT/'candidate_scores.npz',ids=np.array(ids),angles=np.array(ANGLES),candidate_ids=np.array(ids)[candidates],scores=scores)
    sub=np.array([subset(meta[s]) for s in ids]);groups=readcsv(OUT/'input_snapshot/metadata/partition_units.csv')
    gu={g:r['partition_unit_id'] for r in groups for g in r['group_ids'].split(';')};unit_ids=sorted({gu[meta[s]['data_group_id']] for s in ids});uindex={u:i for i,u in enumerate(unit_ids)};units=np.array([uindex[gu[meta[s]['data_group_id']]] for s in ids])
    rng=np.random.default_rng(20261016);weights=rng.multinomial(len(unit_ids),np.full(len(unit_ids),1/len(unit_ids)),size=4000)
    summaries=[]
    for group in ('all','limb','maxilla','mandible','mouse'):
        idx=np.ones(len(ids),bool) if group=='all' else sub==group;size=int(idx.sum());den=np.bincount(units[idx],minlength=len(unit_ids));bootden=weights@den
        for ai,angle in enumerate(ANGLES):
            rank=ranks[idx,ai];positive=scores[idx,ai,0];negative=scores[idx,ai,1:].ravel();n1=len(positive);n0=len(negative);auc=(rankdata(np.r_[positive,negative])[:n1].sum()-n1*(n1+1)/2)/(n1*n0)
            count=np.bincount(units[idx],weights=(rank==1),minlength=len(unit_ids));ci=np.percentile(100*(weights@count)[bootden>0]/bootden[bootden>0],[2.5,97.5])
            summaries.append(dict(subset=group,angle=angle,records=size,top1_count=int((rank==1).sum()),top3_count=int((rank<=3).sum()),top5_count=int((rank<=5).sum()),top1_percent=float(100*(rank==1).mean()),top3_percent=float(100*(rank<=3).mean()),top1_cluster_ci_low=float(ci[0]),top1_cluster_ci_high=float(ci[1]),mean_rank=float(rank.mean()),auc=float(auc)))
    savecsv(OUT/'summary.csv',summaries)
    per=[]
    for qi,sid in enumerate(ids):
        for ai,angle in enumerate(ANGLES):
            slot=1+int(np.argmax(scores[qi,ai,1:]));per.append(dict(sample_id=sid,subset=subset(meta[sid]),source_volume_id=meta[sid]['source_volume_id'],angle=angle,score_matched=float(scores[qi,ai,0]),rank=int(ranks[qi,ai]),best_other_id=ids[candidates[qi,slot]],score_best_other=float(scores[qi,ai,slot]),margin=float(scores[qi,ai,0]-scores[qi,ai,slot])))
    savecsv(OUT/'per_record.csv',per)
    changes=[]
    inputs=readcsv(OUT/'input_file_hashes.csv')
    for r in inputs:
        path=ROOT/r['relative_path'];stat=path.stat()
        if stat.st_size!=int(r['bytes']) or stat.st_mtime_ns!=int(r['mtime_ns']):
            actual=file_sha(path)
            if actual!=r['sha256']:changes.append(dict(path=r['relative_path'],frozen=r['sha256'],current=actual))
    snapshot={n:file_sha(ROOT/n) for n in json.loads((OUT/'protocol.json').read_text())['dataset_fingerprints']}
    savejson(OUT/'input_stability_verification.json',dict(status='passed' if not changes else 'input_changed_after_feature_extraction',input_files=len(inputs),changed_image_files=changes,current_metadata_fingerprints=snapshot,method='Compare size/mtime of all 13,374 inputs; rehash any changed file against the bytes actually scored.'))
    if changes:raise ValueError('Input images changed while running; report frozen results but do not claim current-version completion.')
    savejson(OUT/'verification.json',dict(status='passed',records=len(ids),radiographs=len(ids)*5,image_inputs=len(inputs),candidate_scores=int(scores.size),all_scores_finite=True,all_five_angles_fresh=True,old_scores_reused=False,queries_with_tied_recorded_score=int((scores[:,:,1:]==scores[:,:,0,None]).any(axis=2).sum()),bootstrap='4,000 partition-unit cluster resamples; seed 20261016; conditional on fixed candidate sets and linkage metadata.',analysis_code_sha256=file_sha(__file__)))
    progress('completed',len(ids),len(ids),candidate_scores=int(scores.size),all_five_angles_fresh=True)
    print(json.dumps([s for s in summaries if s['subset']=='all'],indent=2),flush=True)

def main():
    global ROOT,OUT
    ap=argparse.ArgumentParser();ap.add_argument('--dataset',type=Path,default=ROOT);ap.add_argument('--out',type=Path,default=OUT);ap.add_argument('--workers',type=int,default=2);ap.add_argument('--fft-workers',type=int,default=2);args=ap.parse_args()
    ROOT=args.dataset.resolve();OUT=args.out.resolve();assert not OUT.is_relative_to(ROOT)
    protocol,meta,ids,candidates,files,manifest=freeze();prepare_features(meta,ids,files,manifest,args.workers)
    scores=score_all(ids,candidates,args.workers,args.fft_workers);verify_against_reference(ids,candidates,scores);finalize(meta,ids,candidates,scores)
if __name__=='__main__':main()
