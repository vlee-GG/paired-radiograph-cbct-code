"""Render the frozen 320 training/formal composites from released image inputs.

Dependencies: numpy, nibabel, pydicom, Pillow. Original PNG byte reproduction also
requires the same Pillow/font build; selected sections and pixel arrays are fixed.
"""
from pathlib import Path
import csv,argparse,hashlib
import numpy as np
import nibabel as nib
import pydicom
from PIL import Image,ImageDraw,ImageFont

def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig',newline='')))
def main(dataset,out,font_path):
    root=Path(__file__).resolve().parent
    files={r['relative_path']:r for r in read(root/'input_hashes.csv')}
    index={r['sample_id']:r for r in read(dataset/'metadata/metadata.csv')}
    font=ImageFont.truetype(str(font_path),24);out.mkdir(parents=True,exist_ok=True)
    def base(r):
        parent=dataset/'data'/r['species']
        if r['species']=='rabbit':parent=parent/{'maxillary':'maxilla','mandibular':'mandible','limb_unspecified':'limb'}[r['anatomical_region']]
        return parent/r['source_volume_id']/r['sample_id']
    for r in read(root/'display_choices.csv'):
        a=index[r['xray']];b=index[r['crop']]
        xp=base(a)/(a['sample_id']+'__xray__a090.dcm');cp=base(b)/(b['sample_id']+'__cbct_crop.nii.gz')
        for path in [xp,cp]:
            rel=path.relative_to(dataset).as_posix()
            assert hashlib.sha256(path.read_bytes()).hexdigest()==files[rel]['sha256'],rel
        xr=pydicom.dcmread(xp).pixel_array.astype(np.float32)
        cube=np.asanyarray(nib.load(cp).dataobj).astype(np.float32);k=int(r['k_star']);i=int(r['i_star'])
        arrays=[xr,cube.max(axis=2).T,cube[:,:,k].T,cube[i,:,:]]
        captions=['A  90°平片','B  CBCT 最大密度投影','C  CBCT 层面一','D  CBCT 正交层面二']
        canvas=Image.new('L',(1560,1560),0);draw=ImageDraw.Draw(canvas)
        for panel,(arr,caption) in enumerate(zip(arrays,captions)):
            lo,hi=(3.,59745.) if panel==0 else (-1000.0001220703,1780.0089111328)
            aimg=Image.fromarray(np.round(np.clip((arr-lo)/(hi-lo),0,1)*255).astype('uint8'))
            scale=min(700/aimg.width,700/aimg.height);pic=aimg.resize((round(aimg.width*scale),round(aimg.height*scale)),Image.Resampling.LANCZOS)
            x=(panel%2)*780;y=(panel//2)*780;draw.text((x+12,y+10),caption,font=font,fill=240)
            canvas.paste(pic,(x+(780-pic.width)//2,y+55+(700-pic.height)//2))
        target=out/(r['composite_id']+'.png');canvas.save(target,optimize=True)
        print(r['composite_id'], 'original_png_sha256_match',hashlib.sha256(target.read_bytes()).hexdigest()==r['sha256'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dataset',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    p.add_argument('--font',type=Path,default=Path('C:/Windows/Fonts/msyh.ttc'));a=p.parse_args();main(a.dataset,a.out,a.font)
