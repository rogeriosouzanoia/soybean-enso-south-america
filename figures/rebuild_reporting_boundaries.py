from pathlib import Path
from copy import deepcopy
import json,sqlite3,struct,zipfile,re,shutil,difflib,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MPolygon
from shapely.geometry import shape
from shapely.ops import unary_union
from shapely import from_wkb
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
FIG=OUT/'main';FIG.mkdir(exist_ok=True)
SRC=ROOT/'revision'
U=ROOT/'upload'

def polygons(g):
    return [g] if g.geom_type=='Polygon' else list(g.geoms)
def draw(ax,g,fc,ec,lw,z=1):
    for p in polygons(g):
        ax.add_patch(MPolygon(np.asarray(p.exterior.coords),facecolor=fc,edgecolor=ec,lw=lw,zorder=z))

# Delimit reporting groups from complete administrative boundaries, not a hull of soybean units.
states=json.load(open(ROOT/'geometry/ne_50m_admin_1_states_provinces.geojson'))['features']
countries=json.load(open(ROOT/'geometry/ne_50m_admin_0_countries.geojson'))['features']
groups={'Southern Brazil':{'PR','SC','RS'},'MATOPIBA':{'MA','TO','PI','BA'},'Central Brazil':{'AM','DF','GO','MG','MS','MT','PA','RO','RR','SP'}}
col={'Argentina':'#e97864','Southern Brazil':'#2863a8','Central Brazil':'#57aad4','MATOPIBA':'#9356be'}
geoms={name:unary_union([shape(f['geometry']) for f in states if f['properties']['admin']=='Brazil' and f['properties']['postal'] in codes]) for name,codes in groups.items()}
geoms['Argentina']=shape(next(f for f in countries if f['properties']['ADMIN']=='Argentina')['geometry'])
con=sqlite3.connect(SRC/'inputs/retained_administrative_footprint.gpkg')
rows=con.execute('select country,parent_region_name,geom from local_phase_effects group by country,spatial_unit_id').fetchall()
local=[]
for c,p,b in rows:
    env=(b[3]>>1)&7; offset=8+{0:0,1:32,2:48,3:48,4:64}[env]
    local.append((c,p,from_wkb(b[offset:])))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7})
fig=plt.figure(figsize=(2.77,3.06),dpi=200)
ax=fig.add_axes([.015,.015,.97,.97]);ax.set_xlim(-74,-33);ax.set_ylim(-43,7);ax.set_aspect('equal');ax.axis('off')
for f in countries:
    if f['properties']['CONTINENT']=='South America':draw(ax,shape(f['geometry']),'#f5f5f5','#7e7e7e',.35)
for name,g in geoms.items():draw(ax,g,matplotlib.colors.to_rgba(col[name],.13),'none',0,2)
for c,p,g in local:
    name='Argentina' if c=='Argentina' else 'Southern Brazil' if p in ['PR','SC','RS'] else 'MATOPIBA' if p in ['MA','TO','PI','BA'] else 'Central Brazil'
    draw(ax,g,col[name],'none',0,3)
for name,g in geoms.items():draw(ax,g,'none',col[name],.85,4)
labels=[('Argentina',(-64,-35),(-73,-39)),('Southern Brazil',(-52,-28),(-49,-33)),('Central Brazil',(-55,-14),(-73,-10)),('MATOPIBA',(-45,-10),(-47,2))]
for name,xy,tx in labels:
    ax.annotate(name,xy,xytext=tx,color=col[name],fontsize=6.7,fontweight='bold',ha='left',va='center',bbox={'boxstyle':'square,pad=.12','fc':'white','ec':col[name],'lw':.4},arrowprops={'arrowstyle':'-','lw':.5,'color':col[name]},zorder=6)
fig.savefig(FIG/'Reporting_macroregions.pdf');fig.savefig(FIG/'Reporting_macroregions.png',dpi=300);plt.close(fig)
locator=Image.open(FIG/'Reporting_macroregions.png').convert('RGB').resize((277,306),Image.Resampling.LANCZOS)
for n in range(1,6):
    im=Image.open(SRC/f'Figures/Figure_{n}.png').convert('RGB')
    if n==1:im.paste(locator,(1457,1017))
    if n==3:im.paste(locator,(1436,10))
    im.save(FIG/f'Figure_{n}.png',dpi=(300,300))
    im.save(FIG/f'Figure_{n}.tiff',compression='tiff_lzw',dpi=(300,300))

