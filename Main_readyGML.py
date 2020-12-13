#!/usr/bin/env python
# coding: utf-8

# In[4]:


# For data
import networkx as nx
from parse import *
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os
import yaml
import configparser

from itertools import product

# import ipynb.fs.full.SolvePyomo as spyomo
import SolveMIP as mip

def print_answer(status, correct, answ, V, K):
    # print(status)
    if correct == "correct":
        for i, j, k in product(V, V, K):
            if (answ[i][j][k] != None) and (answ[i][j][k] > 0.0):
                print("answer[{}][{}][{}] = {}".format(i, j, k, answ[i][j][k]))
    else:
        print(correct)
    
def set_time(t0, T, Contracts):
    t = {}
    points = set()
    points.add(t0)
    points.add(T)
    for c in Contracts:
        points.add(c['t_s'])
        points.add(c['t_s'] + c['t'])
    points = sorted(list(points))
    for i in range(len(points) - 1):
        e = points[i + 1] - points[i]
        C = set()
        for c in range(len(Contracts)):
            if points[i] >= Contracts[c]['t_s'] and points[i] < Contracts[c]['t_s'] + Contracts[c]['t']:
                C.add(c)
        t[i] = {'em': e, 'Ctm': C}
    return points, t

def unit_measure(strr):
    mmnt = {'b/s':1e-6, 'kb/s': 1e-3, 'Mb/s': 1, 'Gb/s':1e+3, 'Tb/s':1e+6}
    strr = strr.replace(" ", "")
    i = strr.rfind('b/s') + 3
    strr = strr[:i]
    for m in mmnt.keys():
        if m in strr:
            pattern = '{:d}' + m
            k = parse(pattern, strr)
            if k != None:
                return k[0]*mmnt[m]
    return None

random.seed(3)

config = configparser.ConfigParser()
config.read("settings.ini")
TOPO = config["Def"]["TOPO_NAME"]
CONTRACTS_NUM = config["Def"]["ACON_NUM"]
CTYPE1 = config["Def"]["CTYPE1"]
CTYPE2 = config["Def"]["CTYPE2"]
CTYPE3 = config["Def"]["CTYPE3"]
T_MAX = config["Def"]["T_MAX"]
PATH = 'Contracts/{}'.format(TOPO)

H1=nx.read_graphml('topo/{}.graphml'.format(TOPO))
H = nx.Graph(H1)

with open('{}/C={}_TMAX-{}_CT1-{}_CT2-{}_CT3-{}.yaml'.format(PATH, CONTRACTS_NUM, T_MAX, CTYPE1, CTYPE2, CTYPE3), 'r') as f:
	Contracts = yaml.load(f, Loader=yaml.FullLoader)

print("Info about topology:")
print("Name   " + "{}".format(TOPO))
print("Number of nodes " + str(len(H.nodes())))
print("Number of edges " + str(len(H.edges())))

Bandwidth = pd.DataFrame(0, index=np.arange(len(H.nodes())), 
                            columns=np.arange(len(H.nodes())))
BackFlows = pd.DataFrame(0, index=np.arange(len(H.nodes())), 
                            columns=np.arange(len(H.nodes())))

for edge in H.edges(data = True):
    data_edge = H.get_edge_data(edge[0], edge[1])
    bndw = unit_measure(data_edge['LinkLabel'])
    bflw = data_edge['BackgroundFlow']
    Bandwidth[int(edge[1])][int(edge[0])] = bndw
    Bandwidth[int(edge[0])][int(edge[1])] = bndw
    BackFlows[int(edge[1])][int(edge[0])] = bflw
    BackFlows[int(edge[0])][int(edge[1])] = bflw
    

mmin = 10000000
mmax = 0
for c in Contracts:
    if c['t_s'] < mmin:
        mmin = c['t_s']
    if c['t_s'] + c['t'] > mmax:
        mmax = c['t_s'] + c['t']
        
V = set(range(Bandwidth.shape[0])) # число вершин
K = set(range(len(Contracts)))     # количество контрактов К
# t0 = 0                             # начало отсчета t0
# T = 1000                           # конец отсчета в сек
t0 = mmin
T = mmax
timeline, Intervals = set_time(t0, T, Contracts)
M = set(range(len(Intervals)))
print('Number of contracts ')
print(len(Contracts))
print("Timeline len " + str(len(timeline)))
# print(timeline)
print("Time intervals len " + str(len(Intervals)))
# print(Intervals)

print("MIP evaluationg ...")
status, correct_mip, answ_mip, time_mip = mip.evaluate_model(V, K, M, Contracts, Intervals, Bandwidth, BackFlows)
print("Total time: ", time_mip)
print_answer(status, correct_mip, answ_mip, V, K)

# print("Pyomo evaluationg ...")
# status, correct_pyomo, answ_pyomo, time_pyomo = spyomo.evaluate_model(V, K, M, Contracts, Intervals, Bandwidth, BackFlows)
# print("Total time: ", time_pyomo)
# print_answer(status, correct_pyomo, answ_pyomo, V, K, M)

# print("Does MIP and Pyomo answers are equal?")
# print((answ_mip == answ_pyomo).all())


# In[ ]:




