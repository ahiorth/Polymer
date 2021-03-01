### -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 13:50:02 2020

@author: 2902412
"""

#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

Nb_=9-1
Nb_=104
#Nb_=1050-1

def analytical(a,p=Nb_):
    '''
    Model from "Theory of Depolymerization of Long Chain Molecules"
    E. W. Montroll and R. Simha, J. Chem. Phys. 1940
    a = average fraction of bonds cut in molecule
    p = initial number of 
    '''
    Mn=1/(1+a*p)
    Mw=a*a*(1+p)+2*(1-a)*((1-a)**(p+1)+a*(p+1)-1)
    Mw=Mw/(a*a*(1+p)*(1+p))
    DP=Mw/Mn
    return np.array([Mn, Mw,DP])

def ana_distFt(t,a,p=Nb_):
    return np.array(a*t*(1-a)**(t-1)*(2+(p-t)*a)/(p+1))

def ana_distNt(t,a,p=Nb_):
    return np.array(a*(1-a)**(t-1)*(2+(p-t)*a)/(p*a+1))


def plot_Mw(fname='polymer_t.out', xscale='log', yscale='log', xrange=[],yrange=[],x='time'):
    data = pd.read_table(fname,encoding = "ISO-8859-1", sep='\t')
    T=data[x].to_numpy()
    Mn=data['Mn'].to_numpy()
    Mw=data['Mw'].to_numpy()
    PD=data['PD'].to_numpy()
    ana=np.array([])
    if(x=='FractionCut'):
        Tf=np.arange(0,1,0.01)
        ana=analytical(Tf)

    fig, ax1 = plt.subplots()
    if xrange:
        ax1.set_xlim(xrange)
    if yrange:
        ax1.set_ylim(yrange)

    ax2 = ax1.twinx()
    
    #hack
    #Mw=Mw/(1+T)
    l1=ax1.plot(T,Mn/Mn[0],label='Mn')
    l2=ax1.plot(T,Mw/Mw[0],label='Mw')
    lana1=lana2=[]
    if len(ana)>0:
        lana1=ax1.plot(Tf,ana[0],'*',label='Mn-ana')
        lana2=ax1.plot(Tf,ana[1],'^',label='Mw-ana')
        ax2.plot(Tf,ana[2],'x')

    if xscale:
        ax1.set_xscale(xscale)
    if yscale:
        ax1.set_yscale(yscale)
    ax1.set_ylabel('Molecular weight (Normalized)')
    ax2.set_ylabel('Poly Dispersivity')
    ax1.set_xlabel(x)
    l3=ax2.plot(T,PD,c='k',label='PD')
    ax1.grid()
    lns=l1+l2+l3+lana1+lana2

    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs,  bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.show()

def hist_PB(time=1000,fname='hist',ext='.out',p=Nb_,Ft=False):
    fig, axes= plt.subplots(1,1)
    y=np.loadtxt(fname+str(time)+ext)
    if Ft:
        y2=(1+y)*(1+y)
    else:
        y2=y
    n, bins, patches = axes.hist(x=y2, bins=40, color='#0504aa',density=True,
    alpha=0.7, rwidth=0.85)
    print(n)
    print(0.5*(bins[1:]+bins[:-1]))
    #print(patches)
    axes.grid(axis='y', alpha=0.75)
    a=0.0961538
    if not Ft:
        ax2=axes.twinx()
        t=np.arange(0,100,0.1)
        Ft=ana_distFt(t,a)
        Nt=ana_distNt(t,a)
        ax2.plot(t,Ft,c='k')
        ax2.plot(t,Nt,'--',c='r')
        t=0.5*(bins[1:]+bins[:-1])
        ax2.set_ylim([0,max(n)])
        ax2.plot(t,n,'-*')
        ax2.plot(t,t*n/2,'-*')
    t=np.arange(1,100,0.1)
    t=0.5*(bins[1:]+bins[:-1])
    Nt=ana_distNt(t,a)
    #ax2.plot(t,Nt,'--',c='k')
    #ax2.set_yscale([0,:])
    axes.set_xlabel('Number of bounds')
    axes.set_ylabel('Frequency')
 #   axes.text(0.5, 0.5, 'iteration = ' + str(time), horizontalalignment='center',
 #                 verticalalignment='center', transform=axes.transAxes)
    Mn=np.sum(y+1)
    Mw=np.sum((1+y)*(1+y))/Mn
    Mn=Mn/len(y)
    #Mw2=np.mean(57.07*(1+y)*(1+y)/(len(y)))/(1+p)
    Mw2=np.sum((1+y)*(1+y))/np.sum(1+y)/(1+p)
    r=1/p*((p+1)/Mn-1)
    print('alpha=',r, 'Mn=', Mn/(1+p))
    print('Mw modified', Mw2, 'Correct= ', analytical(a))
    print("fraction of bonds at time "+str(time)," ", Mn/(1+p)," ", Mw/(1+p))
    print("average number of bonds at time "+str(time)," ", Mn," ", Mw)

    print("PD=",Mw/Mn)
#    plt.show()
#   plt.close()
#    axes.set_ylim([0,1e-4])


plot_Mw()
plot_Mw(xscale=[],yscale=[],xrange=[1,1000])
plot_Mw(xscale=[],x='FractionCut')
plot_Mw(xscale=[],yscale=[],x='FractionCut')
hist_PB(time=10)
hist_PB(time=10,Ft=True)

y=np.loadtxt('hist10.out')
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

# %%
