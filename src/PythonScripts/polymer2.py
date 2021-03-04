#%%
import matplotlib.pyplot as plt
import numpy as np
from numba import jit

#@jit(nopython=True)
def make_histogram(polymers,nm):
    y=np.append(polymers,nm*[0])
    fig, axes = plt.subplots(1,1)
    n, bins, patches = axes.hist(x=y, bins=30, color='#0504aa',density=True,
                            alpha=0.7, rwidth=0.85)
    axes.grid(axis='y', alpha=0.75)
    axes.set_xlabel('Number of bounds')
    axes.set_ylabel('Frequency')
    plt.show()
#@jit(nopython=True)
def MC_run():    
    p_=104
    N_=10000
    Polymers_=np.array(N_*[p_])
    clock_=0
    time_end_=180000
    dt_=N_
    N_mono_=0
    Tn_=[]
    Mn_=[]
    Mw_=[]
    r_ =[]
    while(clock_<=time_end_ and len(Polymers_)>0):
        if(clock_%dt_==0):
            Tn_.append(clock_)
            Mn_.append(np.sum(Polymers_)/len(Polymers_))
            r_.append(((p_+1)/Mn_[-1]-1)/p_)
            print(Tn_[-1],"\t",Mn_[-1],"\t", Mn_[-1]*57.7,"\t",r_[-1])
            make_histogram(Polymers_,N_mono_)    
        if(len(Polymers_)==1):
            l=0
        else:
            l=np.random.randint(len(Polymers_))
        if(Polymers_[l]>1):
            l1=np.random.randint(Polymers_[l]-1)+1  #from 1...
        else:
            l1=1
        l2=Polymers_[l]-l1
        if(l2==0):
            N_mono_+=1
        else:
            Polymers_=np.append(Polymers_,l2)
        if(l1==1):
            N_mono_ +=1
        else:
            Polymers_=np.append(Polymers_,l1-1)
        Polymers_[l]=0
        Polymers_=Polymers_[Polymers_>0]
        clock_+=1
    return Polymers_,N_mono_

dd,dm=MC_run()

x=np.bincount(np.append(dd,dm*[0]))
np.savetxt('test.out', x, delimiter=',')

# %%
