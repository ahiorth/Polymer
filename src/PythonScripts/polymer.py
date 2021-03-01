#%%
import numpy as np
import matplotlib.pyplot as plt 
import sys
#import numba as nb
#import itertools as it
#from numba import jit
#from numba import int32, float32, int64, boolean    # import the types
#from numba.experimental import jitclass
from celluloid import Camera

ANIM_=False
fig, [axes,axes2] = plt.subplots(1,2)
if ANIM_:     
    
    camera = Camera(fig)
    print('yes')
np.random.seed(2)
#spec = [
#    ('Mw_C', float32),('Mw_H', float32),('Mw_O', float32),('Mw_N', float32),
#    ('Mw_K', float32),('Mw_Na', float32),('Mw_poly', float32),('Mw_mono', float32),
#    ('Nb_mono', int32),('Nb_mono_deg', int32),('Nb_poly', int32),('Nb_polyT', int32[:]),
#    ('N_deg_monomers', int32),('N_deg_step', int32),('fully_degraded', boolean),
#    ('array', float32[:]),          # an array field
#]
#@jitclass(spec)
class polymer:
    def __init__(self,Mw=6e6,Nb=1,Nb_mono_deg=1,NC=3,NH=5,NO=1,NN=1,NK=0,NNa=0):
        """
        solver class for calculating polymer degradation NC, ..., NNa is the number of Carbon,
        hydrogen, oxygen, nitrogen, kalium, and natrium in the polumer monomer
        """
        self.Mw_C    = 12.0107    #g/mol
        self.Mw_H    = 1.0079     #g/mol
        self.Mw_O    = 15.9994    #g/mol
        self.Mw_N    = 28.014     #g/mol
        self.Mw_K    = 39.098     #g/mol
        self.Mw_Na   = 22.989769  #g/mol
        self.Mw_poly = Mw         #g/mol
        self.Mw_mono = NC*self.Mw_C+NH*self.Mw_H+NO*self.Mw_O+NK*self.Mw_K+NNa*self.Mw_Na
        self.Nb_mono = Nb         #number of bindings in monomer connected to the backbone
        self.Nb_mono_deg=Nb_mono_deg
        self.Nb_poly = int(self.Mw_poly/self.Mw_mono*self.Nb_mono)-1
        self.N_polyT=[self.Nb_poly] # a list of the individual length of polymer chain after time T
        self.N_deg_monomers=0     #cummulative number of degraded monomers 
        self.N_deg_step=0         #number of degraded per step

        self.fully_degraded=False

    #@jit(nopython=True)
    def degrade(self):
        """
        pick a random position in a polymer backbone to be degraded
        all elements less or equal to a monomer are removed
        """
        # pick a random degraded polymer,
        
        l=np.random.randint(len(self.N_polyT))
        # ... and place in that polymer
        l1=np.random.randint(self.N_polyT[l]-1)+1  #from 1...
        l2=self.N_polyT[l]-l1
        self.N_polyT[l]=l1-1
        self.N_deg_step=0
        if self.N_polyT[l] <= self.Nb_mono_deg:
            self.N_deg_step+=l1
            self.N_polyT.pop(l) #remove polymer from list
        if l2 <= self.Nb_mono_deg:
            self.N_deg_step+=l
        else:
            self.N_polyT.append(l2) # add to list
        if not self.N_polyT: # list empty, no polymer left
            self.fully_degraded=True
        self.N_deg_monomers += self.N_deg_step
        return self.N_deg_step
    
class polymer_solution:
    def __init__(self, polymer_types, no_polymer_types=[],DT=10000, Tf=1000000):
        """
        takes as input an array of polymer objects, and degrades them
        polymer_types: list of polymer objects
        no_polymer_types: number of each polymer in the solution
        """
        self.no_polymers=self.args_check(polymer_types, no_polymer_types)
        self.polymer_types=polymer_types
        self.active=[i for i in range(len(self.polymer_types))] # all polymers are active
        self.time_end=Tf
        self.DT=DT
        self.clock=0
        self.mass_degraded_t=[]
        self.Na = 6.0221409e+23 # Avogadros number 
        self.Mn=[] # number average molecular weight
        self.Tn=[]


        
    def args_check(self,polymer_types, no_polymer_types):
        if len(no_polymer_types) >0:
            if len(no_polymer_types) == len(polymer_types):
                return no_polymer_types
            elif len(no_polymer_types) ==0:
                return [1]*len(polymer_types)
            else:
                print("The second argument must be empty or have a length equal to number of polymers")
                sys.exit()

    def get_mass_degraded_all(self):
        """
        calculates the mass of degraded polymers in the solution
        """
        mp=0.
        for pol in self.polymer_types:
            mp += pol.N_deg_monomers*pol.Mw_mono
        return mp/self.Na

   # @jit(nopython=True)
    def degrade_polymer(self):
        """
        pick a polymer in solution, degraded it and return the mass of the degraded part
        """
        idx=np.random.randint(len(self.active)) 
        pol=self.polymer_types[self.active[idx]]
        dmp=pol.degrade()
        if pol.fully_degraded:
            self.active.pop(idx)
            print("polymer ", idx, "fully degraded")
        return dmp*pol.Mw_mono/self.Na

    #@jit(nopython=True)
    def degrade_polymer_solution(self):
        while len(self.active)>0 and self.clock<self.time_end:
            if self.clock==0 or self.clock%self.DT ==0:
                self.make_polymer_hist()
            self.mass_degraded_t.append(self.degrade_polymer())
            self.clock +=1
        
        
    def get_polmer_fractions(self):
        Mw=[]
        Mn=0.
        for pol in self.polymer_types:
            if(pol.N_deg_monomers>0):
                dd=[pol.N_deg_monomers]
            else:
                dd=[]
            mm=pol.N_polyT + dd
            Mw += mm #*pol.Mw_mono
            Mn += sum(mm)*pol.Mw_mono/len(mm)
        self.Mn.append(Mn/len(self.polymer_types))
        self.Tn.append(self.clock)
        return Mw
    
    def make_polymer_hist(self):
#        if ANIM_:
        fig, [axes,axes2] = plt.subplots(1,2)
        y=self.get_polmer_fractions()
 #       plt.hist(x=y, bins='auto', color='#0504aa',alpha=0.7, rwidth=0.85)
        n, bins, patches = axes.hist(x=y, bins='auto', color='#0504aa',density=True,
                            alpha=0.7, rwidth=0.85)
        axes.grid(axis='y', alpha=0.75)
        axes.set_xlabel('Number of bounds')
        axes.set_ylabel('Frequency')
        axes.text(0.5, 0.5, 'iteration = ' + str(self.clock), horizontalalignment='center',
                  verticalalignment='center', transform=axes.transAxes)
        #axes.set_ylim([0,1e-4])
        axes2.plot(self.Tn,self.Mn,'-*',c='b')

#        axes2.set_xlim(1000,self.time_end)
#        axes2.set_ylim(1000,max(self.Mn)*10)
        axes2.set_yscale('log')
        axes2.set_xscale('log')
        axes2.grid()
#        maxfreq = n.max()
# Set a clean upper y-axis limit.
#        axes.set_ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
        print(self.clock)
        if ANIM_:
            camera.snap()
        else:
            plt.show()
            plt.close()
            

            
    def plot_results(self):
        t=np.arange(self.clock)   
        cum=np.cumsum(self.mass_degraded_t)
        axes.plot(t,self.mass_degraded_t,label="polymer degraded per time step")
        axes.plot(t,cum,label="cummulative")   

    

        





    

Mw1=6e3 
N1=100   
Tf=1000000
DT=100
#a=list(it.repeat(polymer(Mw=Mw1), N1))
a=[polymer(Mw=Mw1) for i in range(N1)]
sol=polymer_solution(a,DT=DT,Tf=Tf)
sol.degrade_polymer_solution()
#sol.plot_results()
t1=sol.Tn
Mn1=sol.Mn
if(ANIM_):
    animation = camera.animate()
    animation.save('high_molar_weight.gif', writer = 'imagemagick')

#Mw2=1e3
#N2=int(Mw1/Mw2*N1) #start with same mass  
#N2=10000
#print('N2=', N2) 
#a2=[polymer(Mw=Mw2) for i in range(N2)]
#print(len(a2))
#sol2=polymer_solution(a2,DT=DT,Tf=Tf)
#sol2.degrade_polymer_solution()

t2=sol2.Tn
Mn2=sol2.Mn

if(ANIM_):
    animation = camera.animate()
    animation.save('low_molar_weight.gif', writer = 'imagemagick')

plt.close()
plt.plot(t1,Mn1,label=' 6MDa')
plt.plot(t2,Mn2,label=' 1kDa')
plt.yscale('log')
plt.xscale('log')
plt.grid()
plt.legend()

plt.show()
plt.savefig("time_series.png",bbox_inches='tight',transparent=True)

#sol.plot_results()

# %%

# %%
