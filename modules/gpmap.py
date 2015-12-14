from model_simulation import Genome as Genome
from model_simulation import EffectorGene as EffectorGene

def g_p_mapa(host,path):
    dic_tar={}
    for target in host:
        sn=0.0
        for effector in path.eff:
            if target in effector.targets:
                sn+=effector.targets[target]
                if sn>1:
                    sn=1.0
                    break
        dic_tar[target]=sn
    #print dic_tar
    return dic_tar

def g_p_mapa2(host,path,so): #still to do this
    dic_tar={}
    for target in host:
        sn=0.0
        for effector in path:
            if target in path[effector]:
                sn+=path[effector][target]
        dic_tar[target]=sn/(sn+so)
    return dic_tar

def g_p_mapb(host,path):
    for effector in path.eff:
        sn=0.0
        for target in effector.targets:
            if target in host:
                sn+=effector.targets[target]
        effector.g_score=sn
    return 
