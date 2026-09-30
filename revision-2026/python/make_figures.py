"""Create the five descriptive figures for the 2026 TAC revision.

Run python/compute_results.py first, then python python/make_figures.py.
"""
from pathlib import Path
from io import BytesIO
from zipfile import ZipFile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap
ROOT=Path(__file__).resolve().parents[1]
U=ROOT/'data'/'raw'; D=ROOT/'data'/'derived'; O=ROOT/'figures'
O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,
                     'axes.spines.right':False,'savefig.dpi':230,'figure.dpi':140})
COL={'cropland':'#b37335','grassland':'#699b64','forest':'#314d59',
     'broadleaf_clc':'#72975a','coniferous_clc':'#375b51','mixed_clc':'#aa8a50'}
events=pd.read_csv(D/'monthly_event_state_means.csv')
sub=pd.read_csv(D/'corine_event_monthly_state_means.csv')
with ZipFile(U/'august_5km_grid_2003_2018_2022.zip') as z:
    grid={y:pd.read_csv(BytesIO(z.read(f'Germany_monthly_grid_{y}_m8.csv')))
          for y in [2003,2018,2022]}
for y,d in grid.items():
    xy=d.region_id.str.split(',',expand=True).astype(int)
    d['gx']=xy[0]; d['gy']=xy[1]

# Fig 1: actual country footprint sampled by the three stable classes.
d=grid[2018].pivot(index=['gx','gy'],columns='lc_name',values='lc_fraction').fillna(0)
labels=d[['cropland','grassland','forest']].to_numpy().argmax(axis=1)+1
xs=d.index.get_level_values('gx').to_numpy(); ys=d.index.get_level_values('gy').to_numpy()
canvas=np.zeros((ys.max()-ys.min()+1,xs.max()-xs.min()+1),dtype=int)
canvas[ys-ys.min(),xs-xs.min()]=labels
fig,ax=plt.subplots(figsize=(5.4,6.3))
cmap=ListedColormap(['#ffffff',COL['cropland'],COL['grassland'],COL['forest']])
ax.imshow(canvas,origin='lower',interpolation='nearest',cmap=cmap,vmin=0,vmax=3,
          extent=[xs.min()*5,xs.max()*5,ys.min()*5,ys.max()*5])
ax.set_aspect('equal'); ax.set_xlabel('EPSG:3035 easting (km)'); ax.set_ylabel('EPSG:3035 northing (km)')
ax.set_title('Germany study area and sampled 5 km cells',loc='left',weight='bold')
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=COL[c],label=c.capitalize()) for c in ['cropland','grassland','forest']],
          loc='lower left',bbox_to_anchor=(0,-.17),ncol=3,frameon=False)
fig.subplots_adjust(bottom=.14);fig.savefig(O/'Figure_1_AOI.png',bbox_inches='tight');plt.close(fig)

# Fig 2: monthly anomaly trajectories, using original 5-km sampling scale.
fig,axes=plt.subplots(1,3,figsize=(10.8,3.25),sharey=True)
for ax,y in zip(axes,[2003,2018,2022]):
    for c in ['cropland','grassland','forest']:
        a=events[(events.year==y)&(events.lc_name==c)].sort_values('month')
        ax.plot(a.month,a.ndvi_mean,marker='o',ms=3.5,lw=1.8,color=COL[c],label=c.capitalize())
    ax.axhline(0,color='#666666',lw=.7);ax.set_title(str(y),weight='bold')
    ax.set_xticks(range(4,11),['A','M','J','J','A','S','O']);ax.set_xlabel('Month')
axes[0].set_ylabel('State mean NDVI anomaly')
axes[1].legend(ncol=3,frameon=False,loc='lower center',bbox_to_anchor=(.5,-.35))
fig.suptitle('Within-season greenness anomalies in three drought years',x=.05,ha='left',weight='bold')
fig.subplots_adjust(left=.075,right=.99,top=.85,bottom=.25,wspace=.12)
fig.savefig(O/'Figure_2_monthly_events.png',bbox_inches='tight');plt.close(fig)

# Fig 3: July-September month-specific correlations after within-state detrending.
rr=pd.read_csv(D/'monthly_detrended_year_block_correlations.csv')
fig,axes=plt.subplots(1,2,figsize=(8.4,3.45),sharex=True,sharey=True)
for ax,variable,label in zip(axes,['sm_rootzone_anom','vpd_anom_kpa'],['Root-zone soil moisture','Vapour pressure deficit']):
    for c in ['cropland','grassland','forest']:
        a=rr[(rr.predictor==variable)&(rr.land_cover==c)].sort_values('month')
        ax.plot(a.month,a.r,marker='o',ms=3,lw=1.7,color=COL[c],label=c.capitalize())
        if c=='forest': ax.fill_between(a.month.to_numpy(),a.ci_low.to_numpy(),a.ci_high.to_numpy(),color=COL[c],alpha=.15,lw=0)
    ax.axhline(0,color='#777777',lw=.7);ax.set_title(label)
    ax.set_xticks(range(4,11),['A','M','J','J','A','S','O']);ax.set_xlabel('Month')
axes[0].set_ylabel('Pearson r, within state and class detrended')
axes[1].legend(frameon=False,loc='lower left')
fig.suptitle('Calendar-month hydroclimatic associations, 2000–2024',x=.05,ha='left',weight='bold')
fig.subplots_adjust(left=.09,right=.98,top=.82,bottom=.16,wspace=.12)
fig.savefig(O/'Figure_3_monthly_coupling.png',bbox_inches='tight');plt.close(fig)

# Fig 4: CORINE forest subtype event trajectories, separate extraction at 500 m.
fig,axes=plt.subplots(1,3,figsize=(10.8,3.25),sharey=True)
for ax,y in zip(axes,[2003,2018,2022]):
    for c,label in [('broadleaf_clc','Broad-leaved'),('coniferous_clc','Coniferous'),('mixed_clc','Mixed')]:
        a=sub[(sub.year==y)&(sub.lc_name==c)].sort_values('month')
        ax.plot(a.month,a.ndvi_mean,marker='o',ms=3.3,lw=1.8,color=COL[c],label=label)
    ax.axhline(0,color='#777777',lw=.7);ax.set_title(str(y),weight='bold')
    ax.set_xticks(range(4,11),['A','M','J','J','A','S','O']);ax.set_xlabel('Month')
axes[0].set_ylabel('Equal-state mean NDVI anomaly')
axes[1].legend(ncol=3,frameon=False,loc='lower center',bbox_to_anchor=(.5,-.35))
fig.suptitle('Stable CORINE forest classes within MODIS forest',x=.05,ha='left',weight='bold')
fig.subplots_adjust(left=.075,right=.99,top=.85,bottom=.25,wspace=.12)
fig.savefig(O/'Figure_4_forest_subtypes.png',bbox_inches='tight');plt.close(fig)

# Fig 5: maps retain 5-km cell variation that state means suppress.
fig,axes=plt.subplots(1,3,figsize=(11.1,4.2),sharex=True,sharey=True)
for ax,y in zip(axes,[2003,2018,2022]):
    d=grid[y]; f=d[d.lc_name.eq('forest')]
    image=np.full(canvas.shape,np.nan)
    image[f.gy.to_numpy()-ys.min(),f.gx.to_numpy()-xs.min()]=f.ndvi_anom.to_numpy()
    im=ax.imshow(image,origin='lower',interpolation='nearest',cmap='RdBu',norm=TwoSlopeNorm(vmin=-.12,vcenter=0,vmax=.12),
                 extent=[xs.min()*5,xs.max()*5,ys.min()*5,ys.max()*5])
    ax.set_title(str(y),weight='bold');ax.set_xlabel('Easting (km)')
axes[0].set_ylabel('Northing (km)')
cb=fig.colorbar(im,ax=axes,location='bottom',fraction=.055,pad=.12,shrink=.5)
cb.set_label('Forest NDVI anomaly in August')
fig.suptitle('Forest anomaly variation across sampled 5 km cells',x=.05,ha='left',weight='bold')
fig.subplots_adjust(left=.06,right=.99,top=.84,bottom=.22,wspace=.06)
fig.savefig(O/'Figure_5_august_forest_grid.png',bbox_inches='tight');plt.close(fig)
