import itertools as itr
from numpy import random as rk
from model_simulation import Genome as Genome
from model_simulation import EffectorGene as EffectorGene

def NEWHOST(L,K): #length,universe of target
    hx=set()
    while len(hx)<L:
        zx=rk.randint(1,K+1)
        if not zx in hx:
            hx.add(zx)
    hx=list(hx)
    hx.sort()
    return hx

def NEWPATHOGEN(k,neo,nto,t,size): #size here or move it???
#neo is max number of effectors at time 0,nto max numb of targeted genes for each effector
    pth=Genome([],t,size)
    n=rk.randint(1,neo+1) #initial number of effector
    for j in xrange(n):
        lj=rk.randint(1,nto+1) #number of targeted genes
        j=EffectorGene()
        j.targets=dict(itr.izip(NEWHOST(lj,k),rk.random(lj)))
        pth.add_effector(j)
    return pth
#d={effector:{target:score,target:score}}
