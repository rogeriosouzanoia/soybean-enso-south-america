"""Share Figure 1's lower-right cell between i and the exact PCA locator.

Panels a-h are copied from the submitted figure without changing any pixel.
Panel i is laid out from its archived eight means/intervals; no analysis is refit.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image,ImageDraw,ImageFont

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--figure',type=Path,required=True)
    parser.add_argument('--pca',type=Path,required=True)
    parser.add_argument('--cancellation',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(exist_ok=True,parents=True)
    source=Image.open(args.figure).convert('RGB')
    assert source.size==(1734,1598)
    original=np.array(source)
    pca=Image.open(args.pca).convert('RGB')
    assert pca.size==(1724,1523)
    # The map frame and all four coloured labels are copied from Figure 3.
    crop_box=(1436,10,1713,316)
    inset=pca.crop(crop_box)
    inset.save(args.output/'Figure_1j_PCA_inset_source.png')
    data=pd.read_csv(args.cancellation)
    assert len(data)==8
    data['region']=data.macroregion.replace({'MATOPIBA and northeastern frontier':'MATOPIBA'})
    regions=['Argentina','Southern Brazil','Central Brazil','MATOPIBA']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':16,
        'axes.edgecolor':'black','axes.labelcolor':'black','xtick.color':'black',
        'ytick.color':'black','axes.grid':False,'axes.spines.right':False,
        'axes.spines.top':False})
    fig=plt.figure(figsize=(4.10,5.80),dpi=100,facecolor='white')
    ax=fig.add_axes([.31,.25,.63,.62])
    colors={'El Niño':'#cc4f5d','La Niña':'#3574b4'}
    for phase,offset in [('El Niño',.085),('La Niña',-.085)]:
        z=data[data.phase==phase].set_index('region').loc[regions]
        mean=z['mean'].to_numpy()
        ax.errorbar(mean,np.arange(4)+offset,
            xerr=np.array([mean-z.lower.to_numpy(),z.upper.to_numpy()-mean]),
            fmt='o',color=colors[phase],markersize=6.0,linewidth=1.2,
            capsize=0,label=phase)
    ax.set_xlim(0,100);ax.set_xticks([0,25,50,75,100])
    ax.set_ylim(3.45,-.45)
    ax.set_yticks(range(4),['Argentina','Southern\nBrazil','Central\nBrazil','MATOPIBA'])
    ax.tick_params(axis='both',labelsize=15,length=3,pad=3)
    for label in ax.get_yticklabels():label.set_fontweight('bold')
    ax.set_xlabel('Spatial cancellation (%)',fontweight='bold',fontsize=15,labelpad=9)
    fig.text(.015,.97,'i) Spatial cancellation',fontweight='bold',fontsize=18,va='top')
    fig.legend(loc='lower center',bbox_to_anchor=(.50,.04),ncol=2,
        frameon=False,prop={'weight':'bold','size':14},handlelength=1.0,
        handletextpad=.3,columnspacing=.8)
    panel=args.output/'Figure_1i_compact.png'
    fig.savefig(panel,dpi=100);plt.close(fig)
    output=source.copy();draw=ImageDraw.Draw(output)
    # This cell lies entirely to the right of the unchanged maps a-d and below e-h.
    cell=(1048,890,1734,1598)
    draw.rectangle(cell,fill='white')
    output.paste(Image.open(panel).convert('RGB'),(1048,920))
    locator=Image.new('RGB',(277,440),'white');d=ImageDraw.Draw(locator)
    font_path=matplotlib.font_manager.findfont(matplotlib.font_manager.FontProperties(family='DejaVu Sans',weight='bold'))
    font=ImageFont.truetype(font_path,24)
    d.multiline_text((0,0),'j) Reporting\nmacroregions',font=font,spacing=4,fill='black')
    locator.paste(inset,(0,79))
    locator.save(args.output/'Figure_1j_macroregions.png',dpi=(300,300))
    # Original inset pixels are pasted at their original size, without resampling.
    output.paste(locator,(1457,938))
    output.save(args.output/'Figure_1.png',dpi=(300,300))
    after=np.array(output);mask=np.ones(original.shape[:2],dtype=bool)
    mask[cell[1]:cell[3],cell[0]:cell[2]]=False
    assert np.array_equal(original[mask],after[mask])
    assert output.crop((1457,1017,1734,1323)).tobytes()==inset.tobytes()
    report=dict(canvas=list(output.size),unchanged_original_pixels_outside_cell=True,
        original_panels_a_to_h_preserved=True,redrawn_panel_i_values=data.to_dict('records'),
        panel_i_source_sha256=hashlib.sha256(args.cancellation.read_bytes()).hexdigest(),
        locator_source='Figure 3 PCA inset',inset_crop=list(crop_box),
        inset_pixels_preserved_without_resampling=True,
        layout='Panels i and j share the original lower-right cell; no bottom strip')
    (args.output/'Figure_1_layout_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print('Figure 1 compact layout built; a-h and exact inset pixels preserved; i values unchanged.')

if __name__=='__main__':main()
