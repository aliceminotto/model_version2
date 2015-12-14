#!usr/bin/python
import sys
sys.path.append('codes2/')
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
#import modA as mda

class Genome:
    def __init__(self,effectors,time0,size): #effector is a set
        self.eff=effectors
        self.children=set() #changed from dict cause each genome has .t
        self.t=time0
        self.size=size

    def add_effector(self,eff):
        assert eff not in self.eff
        self.eff.add(eff)

    def child(self,genome): #needs new genome
        self.children.add(genome)

    def pop_dyn(self,new_size):
        self.size=new_size

class EffectorGene:
    def __init__(self):
        self.g_score=0.0 #global score
        self.targets={} #dictionary with targets and their score

    def add_target(self,tar,s): #need target and its score, works also to change score
        self.targets[tar]=s

    def remove_target(self,tar): #needs target
        del self.target[tar]
