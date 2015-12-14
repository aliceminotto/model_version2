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
        self.size=[size]
        self.r=None

    def set_r(self,r):
        self.r=r

    def add_effector(self,eff):
        assert eff not in self.eff
        self.eff.append(eff)

    def child(self,genome): #needs new genome
        self.children.add(genome)

    def pop_dyn(self,new_size):
        self.size.append(new_size)

class EffectorGene:
    def __init__(self):
        self.g_score=0.0 #global score
        self.targets={} #dictionary with targets and their score

    def add_target(self,tar,s): #need target and its score, works also to change score
        self.targets[tar]=s

    def remove_target(self,tar): #needs target
        del self.target[tar]

def jumps(j,rngseed,par):
    #par[0] number of jumps
    #par[1] target poool size
    #par[2] max host length
    #par[3] max numb of effector per path at time0
    #par[4] max numb of target per effector gene
    #par[5] mu1,mu2
    #par[6] rates
    T=0 #start of simulation
    strains=set()
    Hn={} #host dictionary
    rng.seed(rngseed)
    for jn in xrange(par[0]):
        l=rng.randint(1,par[2]) #length of host genome
        u=mda.newhost.NEWHOST(l,par[1])
        Hn[jn]=u
    r=0.0
    #j=0
    while r<=.5:
        #j+=1
        pathogen_dic=mda.newhost.NEWPATHOGEN(par[1],par[3],par[4],T,10)
        r=sum(mda.gpmap.g_p_mapa(Hn[0],pathogen_dic).values())/float(len(Hn[0]))
        #print r
        #for effector in pathogen_dic.eff:
        #    print effector.targets.keys()
    pathogen_dic.set_r(r)
    strains.add(pathogen_dic)
    #path_pop={}
    #path_pop[0]=[0,[10]]
    #path_r={}
    #path_r[0]=r
    #path_genomes={}
    #path_genomes[0]=pathogen_dic
    for t in xrange(par[0]*par[7]):
        for el in strains:
            if el.size[-1]>=0:
                print t,el
                el.pop_dyn(mda.population.N_calc(path_pop,path_r,el,path_pop[el][1][-1],par[8])) ####

def main(): #parallelize
    seeds=[]
    NUM_PROCESSES = 1
    SEED=987654321
    children = []
    ##########
    #Parameters
    DT=5000
    NJ=10 #number of jumps
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
