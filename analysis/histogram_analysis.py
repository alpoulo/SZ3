import os
from os import path
import glob

import pandas as pd #make sure installed
import numpy as np #make sure installed
import matplotlib #make sure installed
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Polygon
import seaborn as sns #make sure installed
import math
import sys
import re
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from scipy.stats import pearsonr
from tabulate import tabulate
from functools import reduce

sns.set()
sns.set_style('whitegrid',{'axes.grid':True})
#sns.set_color_codes()
sns.set_context('notebook',font_scale=1.5,rc={'lines.linewidth':2.5,'lines.markeredgewidth':1.,'lines.markersize':10})
#color=['b','g','r','y','m']
font={'family':'serif','weight':'bold','size':14}
matplotlib.rc('font',**font)
#sns.set_palette('tab10')
sns.set(palette='deep')

appinfo={"appname": ["HACC", "CESM2D", "exaalt-copper_dataset1", "ISABEL", "CESM3D", "NYX", "Miranda", "exaalt-copper_dataset2"],
         "dims": [(1073726487), (3600, 1800), (5243, 3137), (500, 500, 100), (3600, 1800, 26), (512, 512, 512), (384, 384, 256), (83, 1077290)],
         "dimsizes": [1, 2, 2, 3, 3, 3, 3, 2]
      }

appdf=pd.DataFrame(data=appinfo)

app="NYX"
errmode="REL"
predmode="lorenzo"
qbins="b16"
dims=tuple(appdf[appdf['appname']==app]['dims'])[0]

if appdf[appdf['appname'] == app]['dimsizes'].item() == 1:
    original_values = dims;
    original_size = original_values * 4;
else:
    original_values = reduce(lambda x, y: x*y, dims)
    original_size = original_values * 4;

bdf=pd.read_pickle('output/' + app + '_baseline_' + errmode + '_' + predmode + '_' + qbins + '_hist.pkl')
adf=pd.read_pickle('output/' + app + '_v5_' + errmode + '_' + predmode + '_' + qbins + '_hist.pkl')

bdf['error_bound'] = pd.to_numeric(bdf['error_bound']);
adf['error_bound'] = pd.to_numeric(adf['error_bound'])

bdf['adaptive_bits'] = pd.to_numeric(bdf['adaptive_bits']);
adf['adaptive_bits'] = pd.to_numeric(adf['adaptive_bits'])

bdf['bits_per_value'] = pd.to_numeric(bdf['compressed_size'] / original_values);
adf['bits_per_value'] = pd.to_numeric(adf['compressed_size'] / original_values);

bdf['perc_quantized'] = (pd.to_numeric(bdf['pred_values'] / original_values)) * 100;
adf['perc_quantized'] = (pd.to_numeric(adf['pred_values'] / original_values)) * 100;
adf['adapt_perc_quantized'] = (pd.to_numeric(adf['adapted_values'] / adf['pred_values'])) * 100;

better_cr = pd.DataFrame(
    columns=['file_name','error_bound','adaptive_bits','bpv_improvement','bpv_improvement(rel%)']
)
better_psnr = pd.DataFrame(
    columns=['file_name','error_bound','adaptive_bits','psnr_improvement','psnr_improvement(rel%)']
)

fnames = bdf['file_name'].unique()
errbounds = bdf['error_bound'].unique()

for fname in fnames:
    tmpbdf = bdf[bdf['file_name'] == fname]
    for eb in errbounds:
        b_bpv = tmpbdf[tmpbdf['error_bound'] == eb]['bits_per_value']
        b_psnr = tmpbdf[tmpbdf['error_bound'] == eb]['psnr']
        tmpadf = adf[(adf['file_name'] == fname) & (adf['error_bound'] == eb)]
        for i in np.arange(1,5):
            diff = (((tmpadf[tmpadf['adaptive_bits'] == i]['bits_per_value']).item()) - b_bpv).item()
            a_bpv = (tmpadf[tmpadf['adaptive_bits'] == i]['bits_per_value']).item()
            if diff < 0:
                better_cr.loc[len(better_cr.index)] = [
                    fname,
                    eb,
                    i,
                    abs(diff),
                    (abs(diff / (b_bpv.item())) * 100)
                ]
            diff = (((tmpadf[tmpadf['adaptive_bits'] == i]['psnr']).item()) - b_psnr).item()
            a_psnr = (tmpadf[tmpadf['adaptive_bits'] == i]['psnr']).item()
            if diff > 0:
                better_psnr.loc[len(better_psnr.index)] = [
                    fname,
                    eb,
                    i,
                    diff,
                    (diff / (b_psnr.item()) * 100)
                ]
#######################################################################################################################################
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.ion()
nbins=64

font = {
#    'family': 'arial',
    'color' : 'black',
    'weight': 'normal',
    'size'  :  8
}

path='plots/' + app + '/' + errmode + '/' + predmode + '/' + qbins;
os.makedirs(path, exist_ok = True)

cr_unique = better_cr.filter(['file_name','error_bound'], axis=1)
cr_unique = cr_unique.drop_duplicates()

ofile = path + '/' + app + '_betterCR' + '.pdf';
with PdfPages(ofile) as pdf:
    for fname,eb in list(zip(cr_unique['file_name'],cr_unique['error_bound'])):
#        print('fname: ' + fname + '\teb: ' + str(eb))
        tmpdf=better_cr[(better_cr['file_name']==fname) & (better_cr['error_bound']==eb)]
        num_better=tmpdf.shape[0]

        bcr=bdf[(bdf['file_name']==fname) & (bdf['error_bound']==eb)]['compression_ratio'].item()
        bpsnr=bdf[(bdf['file_name']==fname) & (bdf['error_bound']==eb)]['psnr'].item()
        bpred=bdf[(bdf['file_name']==fname) & (bdf['error_bound']==eb)]['perc_quantized'].item()
        tmpdict=bdf[(bdf['file_name']==fname) & (bdf['error_bound']==eb)]['quant_histogram'].item()    

        width=(num_better +1) * 4
        fig,axs = plt.subplots(1,num_better+1,sharey=True,figsize=(width,4))
        fig.suptitle(fname + ' ' + str('{:.1e}'.format(eb)))            
        val,weight = zip(*[(k,v) for k,v in tmpdict.items()])     

        axs[0].hist(val,weights=weight,bins=nbins,density=True)
        axs[0].set_title('baseline')
        axs[0].text(0.1, 0.9, 'CR: %.5gx'%bcr, horizontalalignment='left', verticalalignment='center', transform=axs[0].transAxes,fontdict=font)    
        axs[0].text(0.1, 0.85, 'PSNR: %.4g'%bpsnr, horizontalalignment='left', verticalalignment='center', transform=axs[0].transAxes,fontdict=font)        
        axs[0].text(0.7, 0.9, 'SQ: %.2f%%'%bpred, horizontalalignment='left', verticalalignment='center', transform=axs[0].transAxes,fontdict=font)            

        abits=list(better_cr[(better_cr['file_name']==fname) & (better_cr['error_bound']==eb)]['adaptive_bits'])
        for count, ab,ax in zip(np.arange(1,num_better+1),abits,axs.ravel()):
            tmpdict=adf[(adf['file_name']==fname) & (adf['error_bound']==eb) & (adf['adaptive_bits']==ab)]['quant_histogram'].item()
            acr=adf[(adf['file_name']==fname) & (adf['error_bound']==eb) & (adf['adaptive_bits']==ab)]['compression_ratio'].item()  
            apsnr=adf[(adf['file_name']==fname) & (adf['error_bound']==eb) & (adf['adaptive_bits']==ab)]['psnr'].item()  
            apred=adf[(adf['file_name']==fname) & (adf['error_bound']==eb) & (adf['adaptive_bits']==ab)]['perc_quantized'].item()
            aadapt=adf[(adf['file_name']==fname) & (adf['error_bound']==eb) & (adf['adaptive_bits']==ab)]['adapt_perc_quantized'].item()        
            sns.set_theme(style="ticks")
            val,weight = zip(*[(k,v) for k,v in tmpdict.items()])
            axs[count].hist(val,weights=weight,bins=nbins,density=True)
            axs[count].set_title(str(ab) + ' adaptive bits')
            axs[count].text(0.1, 0.9, 'CR: %.5gx'%acr, horizontalalignment='left', verticalalignment='center', transform=axs[count].transAxes,fontdict=font)
            axs[count].text(0.1, 0.85, 'PSNR: %.4g'%apsnr, horizontalalignment='left', verticalalignment='center', transform=axs[count].transAxes,fontdict=font)    
            axs[count].text(0.67, 0.9, 'SQ: %.2f%%'%apred, horizontalalignment='left', verticalalignment='center', transform=axs[count].transAxes,fontdict=font)
            axs[count].text(0.67, 0.85, 'AQ: %.2f%%'%aadapt, horizontalalignment='left', verticalalignment='center', transform=axs[count].transAxes,fontdict=font)            

        pdf.savefig(fig)
        #plt.close()