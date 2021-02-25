### -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 13:50:02 2020

@author: 2902412
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def plot_Mw(fname='polymer_t.out', scale='log', xrange=[],yrange=[]):
    data = pd.read_table(fname,encoding = "ISO-8859-1", sep='\t')
    T=data['time'].to_numpy()
    Mn=data['Mn'].to_numpy()
    Mw=data['Mw'].to_numpy()
    PD=data['PD'].to_numpy()

    fig, ax1 = plt.subplots()
    if xrange:
        ax1.set_xlim(xrange)
    if yrange:
        ax1.set_ylim(yrange)

    ax2 = ax1.twinx()

    l1=ax1.plot(T,Mn,label='Mn')
    l2=ax1.plot(T,Mw,label='Mw')
    if scale:
        ax1.set_xscale(scale)
        ax1.set_yscale(scale)
    ax1.set_ylabel('Molecular weight')
    ax2.set_ylabel('Poly Dispersivity')
    l3=ax2.plot(T,PD,c='k',label='PD')
    ax1.grid()
    lns=l1+l2+l3

    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs,  bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.show()

def hist_PB(time=1000,fname='hist',ext='.out'):
    fig, axes= plt.subplots(1,1)
    y=np.loadtxt(fname+str(time)+ext)
    n, bins, patches = axes.hist(x=y, bins='auto', color='#0504aa',density=True,
                            alpha=0.7, rwidth=0.85)
    axes.grid(axis='y', alpha=0.75)
    axes.set_xlabel('Number of bounds')
    axes.set_ylabel('Frequency')
    axes.text(0.5, 0.5, 'iteration = ' + str(time), horizontalalignment='center',
                  verticalalignment='center', transform=axes.transAxes)
#    axes.set_ylim([0,1e-4])


plot_Mw()
plot_Mw(scale=[],xrange=[1,10000])

hist_PB(time=100000)

#data2.plot(x=' Time (hours)',y='pH-cell',ls=':',ax=ax,grid=True,label='pH-IORCoreSim')
#data3.plot(x='Time',y='pH',ls='--',grid=True,ax=ax,label='pH-1Dsolver')
#data5.plot(x=data.columns[0],y='pH-cell',ls='--',grid=True,ax=ax,label='pH-IORC-NewIP')

#ax2=data2.plot(x=data2.columns[0],y='pH-cell',label='pH-IORCoreSim-chbal')
#data3.plot(x='Time',y='pH',ls='--',grid=True,ax=ax2,label='pH-1Dsolver')
#plt.savefig('mbal_8.png', bbox_inches='tight',transparent=True)
#data4.plot(x='Time',y='pH',ax=ax,grid=True)

#ax2=data3.plot(x='Time',y='Cl-',grid=True)
#data3.plot(x='Time',y=['Na+','Ca+2','Mg+2'],grid=True,ax=ax2)
#ax=data3.plot(x='Time',y=['K+','Cl-','Na+','NO3-','Ca+2'],grid=True)
#ax2=ax.twinx()
#data3.plot(x='Time',y=['pH'],ax=ax2,grid=True)

#p1=data2.plot(x='Time',y=['pH'],label=['pH-charge balance'])
#data3.plot(x='Time',y=['pH'],ax=p1,grid=True)

#ax=data3.plot(x='Time',y='pH',grid=True)
#data4.plot(x='Time',y='pH',ax=ax,grid=True)
