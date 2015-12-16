#!usr/bin/python
import sys
sys.path.append('modules/')
import itertools as itr
from itertools import product
import os,time
from multiprocessing import Process, Queue
import subprocess as sp
from math import floor
import math as mth
import pickle
import numpy as np
from numpy import random as rng
import modA as mda

class Genome:
    def __init__(self,effectors,time0,size): #effector is a set
        self.eff=effectors
        self.children=set() #changed from dict cause each genome has .t
        self.t=time0
        self.size=size
        self.r=None

    def __deepcopy__(self):
        G=Genome([x.__copy__() for x in self.eff],self.t,self.size)
        G.r=self.r
        return G

    def set_r(self,r):
        self.r=r

    def add_effector(self,eff):
        assert eff not in self.eff
        self.eff.append(eff)

    def remove_effector(self,eff):
        assert eff in self.eff
        self.eff.remove(eff)

    def child(self,genome): #needs new genome
        self.children.add(genome)

    def pop_dyn(self,new_size):
        self.size.append(new_size)

class EffectorGene:
    def __init__(self):
        self.g_score=0.0 #global score
        self.targets={} #dictionary with targets and their score

    def __copy__(self):
        E=EffectorGene()
        for target in self.targets:
            E.add_target(target,self.targets[target])
        E.g_score=self.g_score
        return E

    def add_target(self,tar,s): #need target and its score, works also to change score
        self.targets[tar]=s

    def remove_target(self,tar): #needs target
        del self.targets[tar]

def jumps(j,rngseed,par):
    #par[0] number of jumps
    #par[1] target poool size
    #par[2] max host length
    #par[3] max numb of effector per path at time0
    #par[4] max numb of target per effector gene
    #par[5] mu1,mu2
    #par[6] rates
    T=0 #start of simulation
    mda.strains=set()
    Hn={} #host dictionary
    rng.seed(rngseed)
    for jn in xrange(par[0]):
        l=rng.randint(1,par[2]) #length of host genome
        u=mda.newhost.NEWHOST(l,par[1])
        Hn[jn]=u
    r=0.0
    while r<=.5:
        pathogen_dic=mda.newhost.NEWPATHOGEN(par[1],par[3],par[4],T,[10])
        r=sum(mda.gpmap.g_p_mapa(Hn[0],pathogen_dic).values())/float(len(Hn[0]))
    pathogen_dic.set_r(r)
    mda.strains.add(pathogen_dic)
    for jn in xrange(par[0]):
        print 'jump ',jn,' in progress'
        #here change r
        for strain in mda.strains:
            if strain.size[-1]>0:
                r=sum(mda.gpmap.g_p_mapa(Hn[jn],strain).values())/float(len(Hn[jn]))
                #print mda.gpmap.g_p_mapa(Hn[jn],strain).values()
                strain.r=r
            else:
                strain.r=0.0
        print max([el.r for el in mda.strains])
        print Hn[jn]
        #print rng.get_state()
        #for x in mda.strains:
            #for y in x.eff:
                #print y,y.targets.keys(),y.g_score
        state=rng.get_state()
        #print state[-1]
        for t in xrange(par[0]*par[7]/par[0]):
            #print t+(jn*(par[0]*par[7]/par[0]))
            rmax=max([el.r for el in mda.strains])
            flag=0
            toadd=set()
            for el in mda.strains:
                #print mda.strains
                if el.size[-1]>0:
                    new_pth_aux=mda.transformations.transform(el,par[1],par[5],par[4],Hn[jn],par[6],el.size[-1],state)
                    # ^ pathogen,K,[mu1,mu2],nto,host,rates,size
                    raux=sum(mda.gpmap.g_p_mapa(Hn[jn],new_pth_aux).values())/float(len(Hn[jn]))
                    if raux>rmax: #will be add as a new strain
                        toadd.add(new_pth_aux)
                        new_pth_aux.r=raux
                        flag=1
                        el.child(new_pth_aux)

            for el in mda.strains:
                if el.size[-1]>0:
                    el.pop_dyn(mda.population.N_calc(mda.strains,el,el.size[-1],par[8]))

            if flag!=0: #adding new strains to the pool with born time and initial size
                for new_strain in toadd:
                    new_strain.t=t+(jn*(par[0]*par[7]/par[0]))
                    new_strain.size=[10]
                    mda.strains.add(new_strain)

        print max([el.size[-1] for el in mda.strains])
        pickle.dump([[el.t,el.size] for el in mda.strains], open("testpops"+str(jn)+".p", "wb"),protocol=2)
        pickle.dump([el.r for el in mda.strains], open("r"+str(jn)+".p", "wb"),protocol=2)

    print("test completed")

def main(): #parallelize
    seeds=[]
    NUM_PROCESSES = 1
    SEED=987654321
    children = []
    ##########
    #Parameters
    DT=50
    NJ=5 #number of jumps
    K=100 #target pool size
    LHmax=50 #max host genome length
    NEO=5 #max numer of effector per pathogen at time 0
    NTO=10 #max number of target for each effector gene
    MU1=1.0/3
    MU2=2.0/3
    m1=0.01 #mutation
    m2=0.01 #duplication
    m3=0.01 #deletion
    m4=0.01 #hgt
    NH=10**5 #host population
    PV=[NJ,K,LHmax,NEO,NTO,[MU1,MU2],[m1,m2,m3,m4],DT,NH] #parameters values
    ##########
    for process in xrange(NUM_PROCESSES):
        pid = os.fork()
        seeds.append(SEED-process)
        if pid:
            children.append(pid)
        else:
            i=str(process)+".dat"
            jumps(i,seeds[process],PV)
            os._exit(0)

    for i, child in enumerate(children):
        os.waitpid(child, 0)

if __name__ == "__main__":
    main()
