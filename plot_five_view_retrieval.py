"""Figure 6: cumulative ranks and anatomical-subset Top-1, all release records."""
from pathlib import Path
import argparse,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dataset',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--scores',type=Path);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    z=np.load(args.scores or args.dataset/'validation/multiview_retrieval/candidate_scores.npz',allow_pickle=False);scores=z['scores'];ids=z['ids'];angles=z['angles'];rank=1+(scores[:,:,1:]>scores[:,:,0,None]).sum(axis=2)
    meta={r['sample_id']:r for r in csv.DictReader((args.dataset/'metadata/metadata.csv').open(encoding='utf-8-sig'))}
    subset=lambda s:'Mouse' if meta[s]['species']=='mouse' else {'limb_unspecified':'Rabbit limb','maxillary':'Rabbit maxilla','mandibular':'Rabbit mandible'}[meta[s]['anatomical_region']]
    groups=['Rabbit limb','Rabbit maxilla','Rabbit mandible','Mouse'];heat=np.array([[(rank[np.array([subset(s)==g for s in ids]),ai]==1).mean()*100 for ai in range(5)] for g in groups]);nn=[sum(subset(s)==g for s in ids) for g in groups]
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'axes.labelsize':7,'axes.titlesize':7,'legend.fontsize':6.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6})
    fig,(ax,bx)=plt.subplots(1,2,figsize=(183/25.4,79/25.4),gridspec_kw={'width_ratios':[1.15,1]});fig.subplots_adjust(left=.075,right=.90,bottom=.19,top=.87,wspace=.54)
    colors=['#0072B2','#56B4E9','#202020','#D55E00','#CC79A7'];styles=['-','--','-',':','-.'];x=np.arange(1,21)
    for ai,a in enumerate(angles):ax.plot(x,[(rank[:,ai]<=k).mean()*100 for k in x],label=f'{a}°',color=colors[ai],ls=styles[ai],lw=1.25)
    ax.plot(x,x*5,color='.65',lw=.8,ls='--',label='Chance');ax.set(xlim=(1,20),ylim=(0,102),xlabel='Candidate rank, k',ylabel='Recorded assignment within top k (%)',xticks=[1,3,5,10,15,20],yticks=[0,25,50,75,100]);ax.legend(loc='lower right',ncol=2,columnspacing=.9,handlelength=2)
    im=bx.imshow(heat,vmin=0,vmax=100,cmap='cividis',aspect='auto');bx.set(xticks=range(5),xticklabels=[f'{a}°' for a in angles],yticks=range(4),yticklabels=[f'{g}\n(n = {n:,})' for g,n in zip(groups,nn)],xlabel='Nominal radiographic angle');bx.tick_params(length=0,pad=5)
    for j in range(4):
        for i in range(5):bx.text(i,j,f'{heat[j,i]:.1f}',ha='center',va='center',fontsize=6.6,color='white' if heat[j,i]<45 else 'black')
    cb=fig.colorbar(im,ax=bx,fraction=.045,pad=.04);cb.set_label('Top-1 retrieval (%)',fontsize=7);cb.ax.tick_params(labelsize=6);cb.set_ticks([0,50,100])
    ax.text(-.17,1.08,'a',transform=ax.transAxes,fontweight='bold',fontsize=9);bx.text(-.65,1.08,'b',transform=bx.transAxes,fontweight='bold',fontsize=9)
    for ext in ('pdf','svg','png'):fig.savefig(args.out/f'Figure_6.{ext}',dpi=300)
    with (args.out/'Figure_6_source_data.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['panel','subset','nominal_angle_deg','k','percent','records'])
        for ai,a in enumerate(angles):
            for k in x:w.writerow(['a','all',a,k,100*(rank[:,ai]<=k).mean(),len(ids)])
        for j,g in enumerate(groups):
            for i,a in enumerate(angles):w.writerow(['b',g,a,1,heat[j,i],nn[j]])
    print('Figure 6 exported; records',len(ids),'angles',angles.tolist())
if __name__=='__main__':main()

