import itertools as itr
from numpy import random as rk
import modA as mda
import math
from model_simulation import Genome as Genome
from model_simulation import EffectorGene as EffectorGene

#muation,deletion,dulication,hgt,add-lose target
'''def deepcopy(dictionary):
    new_dic={}
    for key in dictionary:
        new_dic[key]=dictionary[key].copy()
    return new_dic'''

def tgain(pathogen,eff,K):
    #changing pathogen
    zk=rk.randint(1,K+1)
    while zk in eff.targets:
        zk=rk.randint(1,K+1)
    eff.add_target(zk,rk.random())
    return pathogen

def tremove(pathogen,eff):
    #changing pathogen
    removethis=rk.choice(eff.targets.keys())
    eff.remove_target(removethis)
    return pathogen

def mutation(pathogen,eff,K,mu):
    #new_pathogen=pathogen.__deepcopy__() #we change new_pathogen from now on
    #print mda.strains #needs to be update accordingly
    #ind=pathogen.eff.index(eff)
    #new_eff=new_pathogen.eff[ind]
    for target in eff.targets:
        y=rk.randn()
        eff.targets[target]+=y
        if eff.targets[target]<0.0:
            eff.targets[target]=0.0
    ranx=rk.random()
    if ranx<mu[0]:
        pathogen=tgain(pathogen,eff,K)
    elif ranx>=mu[1]:
        pathogen=tremove(pathogen,eff)
    return pathogen

def deletion(pathogen,eff):
    #new_pathogen=pathogen.__deepcopy__() #twe change new_pathogen from now on
    #ind=pathogen.eff.index(eff)
    #new_eff=new_pathogen.eff[ind]
    pathogen.remove_effector(eff)
    return pathogen

def duplication(pathogen,eff):
    #new_pathogen=pathogen.__deepcopy__() #we change new_pathogen from now on
    #ind=pathogen.eff.index(eff)
    #new_eff=new_pathogen.eff[ind]
    pathogen.add_effector(eff.__copy__())
    return pathogen

def hgt(pathogen,nto,k):
    #new_pathogen=pathogen.__deepcopy__() #we change new_pathogen from now on
    new_eff=EffectorGene()
    lj=rk.randint(1,nto+1)  #number of targeted genes
    new_eff.targets=dict(itr.izip(mda.newhost.NEWHOST(lj,k),rk.random(lj)))
    pathogen.add_effector(new_eff)
    return pathogen

def probabilities(pathogen,host,rates,pop):
    dic_prob={} #keys are eff, values are probabilities
    #print pathogen
    #print mda.gpmap.g_p_mapa(host,pathogen)
    av_score_tar=sum(mda.gpmap.g_p_mapa(host,pathogen).values())/float(len(host))
    #print av_score_tar
    mda.gpmap.g_p_mapb(host,pathogen)
    for effector in pathogen.eff:
        dic_prob[effector]=[rates[0]*pop,rates[1]*(1-av_score_tar)*pop,
                            (rates[2]*pop)*math.exp(-effector.g_score),rates[3]*pop/len(pathogen.eff),
                            pop*((1-rates[0])+(1-rates[1])*(1-av_score_tar)+1-(rates[2]*
                            math.exp(-effector.g_score))+1-(rates[3]/len(pathogen.eff)))]
    return dic_prob

def events(pathogen,host,rates,pop):
    dic_prob=probabilities(pathogen,host,rates,pop)
    appening={}
    for effector in dic_prob:
        k=rk.random()
        sum_prob=sum(dic_prob[effector])
        n=0
        alfa=k*sum_prob
        som=0.0
        for i in dic_prob[effector]:
            #print effector
            som+=i
            if som>=alfa:
                break
            n+=1
        appening[effector]=n
    return appening

def transform(pathogen,K,mu,nto,host,rates,pop,state):
    rk.set_state(state)
    appening=events(pathogen,host,rates,pop)
    #print len(pathogen.eff)
    #print appening.values()
    effectors_past=list(pathogen.eff)
    new_pathogen=pathogen.__deepcopy__() #we change new_pathogen from now on
    effector_to_change=list(new_pathogen.eff)
    for effector in effector_to_change: #they are list and are ordered
        ind=effector_to_change.index(effector)
        eff_old=pathogen.eff[ind]
        event=appening[eff_old]
        if event==4:
            pass
        elif event==0:
            #print 0
            new_pathogen=mutation(new_pathogen,effector,K,mu)
        elif event==1:
            #print 1
            new_pathogen=duplication(new_pathogen,effector)
        elif event==2:
            #print 2
            new_pathogen=deletion(new_pathogen,effector)
        elif event==3:
            #print 3
            new_pathogen=hgt(new_pathogen,nto,K)
    return new_pathogen
