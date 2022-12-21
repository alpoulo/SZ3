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



app="CESM2D"
errmode="ABS"
predmode="lorenzo"
qbins="b16"
#dims=[3600,1800,26]
dims=[512,512,512]

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
#adf['bits_per_value'] = pd.to_numeric(adf['bits_per_value'])

better_cr = pd.DataFrame(
    columns=['file_name','error_bound','adaptive_bits','psnr','bits_per_value','bpv_improvement','bpv_improvement(rel%)']
)
better_psnr = pd.DataFrame(
    columns=['file_name','error_bound','adaptive_bits','psnr','bits_per_value','psnr_improvement','psnr_improvement(rel%)']
)

fnames = bdf['file_name'].unique()
errbounds = bdf['error_bound'].unique()

for fname in fnames:
    tmpbdf = bdf[bdf['file_name'] == fname]
    for eb in errbounds:
        better_cr_flag = 0;
        better_psnr_flag = 0;
        b_bpv = tmpbdf[tmpbdf['error_bound'] == eb]['bits_per_value']
        b_psnr = tmpbdf[tmpbdf['error_bound'] == eb]['psnr']
        tmpadf = adf[(adf['file_name'] == fname) & (adf['error_bound'] == eb)]
        for i in np.arange(1,5):
            diff = (((tmpadf[tmpadf['adaptive_bits'] == i]['bits_per_value']).item()) - b_bpv).item()
            a_bpv = (tmpadf[tmpadf['adaptive_bits'] == i]['bits_per_value']).item()
            a_psnr = (tmpadf[tmpadf['adaptive_bits'] == i]['psnr']).item()
            if diff < 0:
                better_cr.loc[len(better_cr.index)] = [
                    fname,
                    eb,
                    i,
                    a_psnr,
                    a_bpv,
                    abs(diff),
                    (abs(diff / (b_bpv.item())) * 100)
                ]
                better_cr_flag = 1;
            diff = (((tmpadf[tmpadf['adaptive_bits'] == i]['psnr']).item()) - b_psnr).item()
            a_psnr = (tmpadf[tmpadf['adaptive_bits'] == i]['psnr']).item()
            if diff > 0:
                better_psnr.loc[len(better_psnr.index)] = [
                    fname,
                    eb,
                    i,
                    a_psnr,
                    a_bpv,
                    diff,
                    (diff / (b_psnr.item()) * 100)
                ]
                better_psnr_flag = 1;
        if (better_psnr_flag):
            better_psnr.loc[len(better_psnr.index)] = [
                fname,
                eb,
                0,
                b_psnr.item(),
                b_bpv.item(),
                0,
                0
            ]            
        if (better_cr_flag):
            better_cr.loc[len(better_cr.index)] = [
                fname,
                eb,
                0,
                b_psnr.item(),
                b_bpv.item(),
                0,
                0
            ]   

better_cr=better_cr.sort_values(by=['file_name','error_bound','adaptive_bits'],ascending=[True,False,True])
better_psnr=better_psnr.sort_values(by=['file_name','error_bound','adaptive_bits'],ascending=[True,False,True])