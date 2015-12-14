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
import modA as mda

class Genome:
    def __init__(self,effectors):
        self.eff=effectors
        self.eff=[]
        self.children={} #child:time

    def add_effector(eff):
        self.eff.append(eff)

    def child(genome,t): #needs new genome and time of creation
        self.children[genome]=t

class EffectorGene:
    def __init__(self):
        self.g_score=0.0 #global score
        self.targets={} #dictionary with targets and their score

    def add_target(tar,s): #need target and its score, works also to change score
        self.targets[tar]=s

    def remove_target(tar): #needs target
        del self.target[tar]
