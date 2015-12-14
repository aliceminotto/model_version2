import itertools as itr
from numpy import random as rk

#muation,deletion,dulication,hgt,add-lose target
def deepcopy(dictionary):
    new_dic={}
    for key in dictionary:
        new_dic[key]=dictionary[key].copy()
    return new_dic

def tgain(pathogen,eff,K):
    new_pathogen=deepcopy(pathogen)
    zk=rk.randint(1,K+1)
    while zk in new_pathogen[eff]:
        zk=rk.randint(1,K+1)
    new_pathogen[eff][zk]=rk.random()
    return new_pathogen

def tremove(pathogen,eff):
    new_pathogen=deepcopy(pathogen)
    removethis=rk.choice(new_pathogen[eff].keys())
    del new_pathogen[eff][removethis]
    return new_pathogen

def mutation(pathogen,eff,K,mu1,mu2):
    new_pathogen=deepcopy(pathogen)
    for target in new_pathogen[eff]:
        y=rk.randn()
        new_pathogen[eff][target]+=y
        if new_pathogen[eff][target]<0.0:
            new_pathogen[eff][target]=0.0
    ranx=rk.random()
    if ranx<mu1:
        new_pathogen=tgain(new_pathogen,eff,K)
    elif ranx>=mu2:
        new_pathogen=tremove(new_pathogen,eff)
    return new_pathogen

def deletion(pathogen,eff):
    new_pathogen={}
    for key in pathogen:
        if key!=eff:
            new_pathogen[key]=pathogen[key].copy()
    return new_pathogen

def duplication(pathogen,eff):
    new_pathogen=deepcopy(pathogen)
    l=max(new_pathogen.keys())+1
    new_pathogen[l]=pathogen[eff].copy()
    return new_pathogen

def hgt(pathogen,newhost,nto,k):
    new_pathogen=deepcopy(pathogen)
    l=max(new_pathogen.keys())+1
    new_pathogen[l]={}
    lj=rk.randint(1,nto+1) #number of targeted genes
    new_pathogen[l]=dict(itr.izip(newhost.NEWHOST(lj,k),rk.random(lj)))
    return new_pathogen
