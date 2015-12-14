from model_simulation import Genome as Genome
from model_simulation import EffectorGene as EffectorGene

def N_calc(strains,el,N_old,NH):
    sommatoria=0
    Ni=N_old
    sommatoria=sum([x.size[-1]*x.r for x in strains])
    print sommatoria
    Ni_t1=((el.r*NH*Ni)/((0.001*NH)+sommatoria))#-Ni
    if Ni_t1<1:
        Ni_t1=0
    return Ni_t1
