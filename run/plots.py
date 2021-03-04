### -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 13:50:02 2020

@author: 2902412
"""

#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

Nb_=9-1
Nb_=104
a_=0.0961538
Time_=1000

#Nb_=105131
#a_=0.0951194
#Time_=1000000
try:
    os.chdir(str(Nb_))
except:
    print('Could not change directory')

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
        Tf=np.arange(0,1,0.001)
        ana=analytical(Tf)

    fig, ax1 = plt.subplots()
    if xrange:
        ax1.set_xlim(xrange)
    if yrange:
        ax1.set_ylim(yrange)

    ax2 = ax1.twinx()
    
    #hack
    #Mw=Mw/(1+T)
    l1=ax1.plot(T,Mn/Mn[0],'*',c='b',label='Mn')
    l2=ax1.plot(T,Mw/Mw[0],'^',c='r',label='Mw')
    lana1=lana2=[]
    if len(ana)>0:
        lana1=ax1.plot(Tf,ana[0],'-', c='b')
        lana2=ax1.plot(Tf,ana[1],'-', c='r')
        ax2.plot(Tf,ana[2],'-',c='g')

    if xscale:
        ax1.set_xscale(xscale)
    if yscale:
        ax1.set_yscale(yscale)
    ax1.set_ylabel('Molecular weight (Normalized)')
    ax2.set_ylabel('Polydispersity')
    ax1.set_xlabel(x)
    l3=ax2.plot(T,PD,'x',c='g',label='PDI')
    ax1.grid()
    lns=l1+l2+l3+lana1+lana2

    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs,  bbox_to_anchor=(1.08, 1), loc='upper left')
    if yscale:
        fname='../../fig/polymer_'+str(Nb_)+'.png'
        plt.savefig(fname, bbox_inches='tight',transparent=True)
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
    a=a_
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

def hist_PB2(time=1000,fname='denst',ext='.out',p=Nb_):
    fig, axes= plt.subplots(1,1)
    Ndenst = np.loadtxt(fname+str(time)+ext)
    M=np.array([i+1 for i in range(len(Ndenst))])
    Mdenst=M*Ndenst
    axes.plot(M,Ndenst/np.sum(Ndenst), '^', c='b')
    axes.plot(M,Mdenst/np.sum(Mdenst),  '*',c='r')
    
    a=a_
    t=np.arange(0,100,0.001)
    Ft=ana_distFt(t,a)
    Nt=ana_distNt(t,a)
    axes.plot(t,Ft,c='r', label='$F_t$')
    axes.plot(t,Nt,c='b', label='$N_t$')
    axes.set_xlabel('Number of chains')
    axes.set_ylabel('Frequency')
    axes.set_xscale('log')
    axes.grid()
    axes.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    fname='../../fig/polymer_hist'+str(Nb_)+'.png'
    plt.savefig(fname, bbox_inches='tight',transparent=True)

 #   axes.text(0.5, 0.5, 'iteration = ' + str(time), horizontalalignment='center',
 #                 verticalalignment='center', transform=axes.transAxes)


def plotMwt(Mwa=[6e3],M0=57.07,k=1,norm=False, savefig=False, data=False):
    fig, ax1 = plt.subplots()
    ln=[]
    lns=[]
    ax2 = ax1.twinx()
    for idx,Mw0 in enumerate(Mwa):
        ax1.set_xlabel('Time [kt]')
        p = int(Mw0/M0-1)
        if(data):
            t=np.logspace(-1, 4, num=1000)
        else:
            t=np.logspace(-7, 4, num=1000)

        a=1-np.exp(-k*t)
        Mn,Mw,DPI=analytical(a,p=p)
        if not norm:
            Mn *= M0*(1+p)
            Mw *= M0*(1+p)
        if(idx==0):
            ln.append(ax1.plot(t,Mn,c='b',label='Mn'))
            ln.append(ax1.plot(t,Mw,c='r',label='Mw'))
            ln.append(ax2.plot(t,DPI,c='k',label='PDI'))
        else:
            ax1.plot(t,Mn,c='b')
            ax1.plot(t,Mw,c='r')
            ax2.plot(t,DPI,c='k')

            
    if data:
        folder='/mnt/c/cygwin64/home/2902412/GitHub/Polymer/data/'
        datafiles=['pam10kDa.txt','pam600KDa.txt','pam6MDa.txt']
        label=['10kDa','600kDa','6MDa']
        Tdata=[]
        Mwdata=[]
        symbol=['s', 'D', 'o']
        col=['b','g','r']
        for idx,files in enumerate(datafiles):
            fname = folder + files
            df = pd.read_table(fname,encoding = "ISO-8859-1", sep='\t')
            Tdata.append(df['Time'].to_numpy())
            Mwdata.append(df['Mw'].to_numpy())
            ln.append(ax1.plot(Tdata[-1],Mwdata[-1],symbol[idx],label=label[idx],
            markersize=10,c=col[idx]))
        ax1.set_xlabel('Time [Days]')

    for ll in ln:
            lns += ll
    labs = [l.get_label() for l in lns]
    print(lns)
    ax1.legend(lns, labs,  bbox_to_anchor=(1.1, 1), loc='upper left')
    ax1.set_ylabel('Molecular weight')
    ax2.set_ylabel('Polydispersity')
   
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.grid()
    if savefig:
        if(data):
            fname='../../fig/polymer_t'+str(p)+'data.png'
        else:    
            fname='../../fig/polymer_t'+str(p)+'.png'
        plt.savefig(fname, bbox_inches='tight',transparent=True)
    plt.show()
    print(p)
#plot_Mw(xscale=[],x='FractionCut')
#plot_Mw(xscale=[],yscale=[],x='FractionCut')
#hist_PB2(time=Time_)

#plotMwt()
plotMwt(Mwa=[8039,800e3,3e6])
plotMwt(Mwa=[8039,800e3,3e6],M0=72.5, data=True,k=800e-7,savefig=True)

#data = pd.read_table('polymer_t.out',encoding = "ISO-8859-1", sep='\t')
#T=data['time'].to_numpy()
#alpha=data['FractionCut'].to_numpy()

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
