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
from itertools import product

import ipynb.fs.full.SolvePyomo as spyomo
import ipynb.fs.full.SolveMIP as mip

def print_answer(status, correct, answ, V, K, M):
    print(status)
    if correct == "correct":
        for i, j, k, t in product(V, V, K, M):
            if (answ[i][j][k][t] != None) and (answ[i][j][k][t] > 0.0):
                print("answer[{}][{}][{}][{}] = {}".format(i, j, k, t, answ[i][j][k][t]))
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

H1=nx.read_graphml('change_topo/cBren.graphml')
H = nx.Graph(H1)

print("Info about topology:")
print("Name   " + "cBren")
print("Number of nodes" + str(H.nodes()))
print("Number of edges" + str(H.edges()))

Bandwidth = pd.DataFrame(0, index=np.arange(len(H.nodes())), 
                            columns=np.arange(len(H.nodes())))
BackFlows = pd.DataFrame(0, index=np.arange(len(H.nodes())), 
                            columns=np.arange(len(H.nodes())))

for edge in H.edges(data = True):
    data_edge = H.get_edge_data(edge[0], edge[1])
    bndw = unit_measure(data_edge['LinkLabel'])
    bflw = data_edge['BackgroundFlow']
#     print(Bandwidth.shape)
#     print(edge)
#     print(Bandwidth[0][237])
    Bandwidth[int(edge[1])][int(edge[0])] = bndw
    Bandwidth[int(edge[0])][int(edge[1])] = bndw
    BackFlows[int(edge[1])][int(edge[0])] = bflw
    BackFlows[int(edge[0])][int(edge[1])] = bflw
    
    
Contracts = [{'volume': 10, 't_s': 0, 't': 10, 'A': 16, 'B': 36},
 {'volume': 10, 't_s': 5, 't': 10, 'A': 18, 'B': 36},
 {'volume': 10, 't_s': 10, 't': 10, 'A': 9, 'B': 36},
 {'volume': 10, 't_s': 15, 't': 10, 'A': 10, 'B': 36},
 {'volume': 10, 't_s': 20, 't': 10, 'A': 12, 'B': 36},
 {'volume': 10, 't_s': 25, 't': 10, 'A': 13, 'B': 36},
 {'volume': 10, 't_s': 30, 't': 10, 'A': 16, 'B': 27},
 {'volume': 10, 't_s': 35, 't': 10, 'A': 16, 'B': 31},
 {'volume': 10, 't_s': 40, 't': 10, 'A': 16, 'B': 32},
 {'volume': 10, 't_s': 45, 't': 10, 'A': 16, 'B': 33},
 {'volume': 10, 't_s': 50, 't': 10, 'A': 16, 'B': 35},
 {'volume': 10, 't_s': 55, 't': 10, 'A': 17, 'B': 36},
 {'volume': 10, 't_s': 60, 't': 10, 'A': 18, 'B': 27},
 {'volume': 10, 't_s': 65, 't': 10, 'A': 18, 'B': 31},
 {'volume': 10, 't_s': 70, 't': 10, 'A': 18, 'B': 32},
 {'volume': 10, 't_s': 75, 't': 10, 'A': 18, 'B': 33},
 {'volume': 10, 't_s': 80, 't': 10, 'A': 18, 'B': 35},
 {'volume': 10, 't_s': 85, 't': 10, 'A': 23, 'B': 36},
 {'volume': 10, 't_s': 90, 't': 10, 'A': 24, 'B': 36},
 {'volume': 10, 't_s': 95, 't': 10, 'A': 26, 'B': 36},
 {'volume': 10, 't_s': 100, 't': 10, 'A': 4, 'B': 16},
 {'volume': 10, 't_s': 105, 't': 10, 'A': 4, 'B': 18},
 {'volume': 10, 't_s': 110, 't': 10, 'A': 4, 'B': 36},
 {'volume': 10, 't_s': 115, 't': 10, 'A': 5, 'B': 16},
 {'volume': 10, 't_s': 120, 't': 10, 'A': 5, 'B': 18},
 {'volume': 10, 't_s': 125, 't': 10, 'A': 5, 'B': 36},
 {'volume': 10, 't_s': 130, 't': 10, 'A': 6, 'B': 16},
 {'volume': 10, 't_s': 135, 't': 10, 'A': 6, 'B': 18},
 {'volume': 10, 't_s': 140, 't': 10, 'A': 6, 'B': 36},
 {'volume': 10, 't_s': 145, 't': 10, 'A': 7, 'B': 16},
 {'volume': 10, 't_s': 150, 't': 10, 'A': 7, 'B': 18},
 {'volume': 10, 't_s': 155, 't': 10, 'A': 7, 'B': 36},
 {'volume': 10, 't_s': 160, 't': 10, 'A': 9, 'B': 12},
 {'volume': 10, 't_s': 165, 't': 10, 'A': 9, 'B': 13},
 {'volume': 10, 't_s': 170, 't': 10, 'A': 9, 'B': 16},
 {'volume': 10, 't_s': 175, 't': 10, 'A': 9, 'B': 27},
 {'volume': 10, 't_s': 180, 't': 10, 'A': 9, 'B': 31},
 {'volume': 10, 't_s': 185, 't': 10, 'A': 9, 'B': 32},
 {'volume': 10, 't_s': 190, 't': 10, 'A': 9, 'B': 33},
 {'volume': 10, 't_s': 195, 't': 10, 'A': 9, 'B': 35},
 {'volume': 10, 't_s': 200, 't': 10, 'A': 10, 'B': 27},
 {'volume': 10, 't_s': 205, 't': 10, 'A': 10, 'B': 31},
 {'volume': 10, 't_s': 210, 't': 10, 'A': 10, 'B': 32},
 {'volume': 10, 't_s': 215, 't': 10, 'A': 10, 'B': 33},
 {'volume': 10, 't_s': 220, 't': 10, 'A': 10, 'B': 35},
 {'volume': 10, 't_s': 225, 't': 10, 'A': 11, 'B': 36},
 {'volume': 10, 't_s': 230, 't': 10, 'A': 12, 'B': 18},
 {'volume': 10, 't_s': 235, 't': 10, 'A': 12, 'B': 27},
 {'volume': 10, 't_s': 240, 't': 10, 'A': 12, 'B': 31},
 {'volume': 10, 't_s': 245, 't': 10, 'A': 12, 'B': 32},
 {'volume': 10, 't_s': 250, 't': 10, 'A': 12, 'B': 33}]
 # {'volume': 10, 't_s': 255, 't': 10, 'A': 12, 'B': 35},
 # {'volume': 10, 't_s': 260, 't': 10, 'A': 13, 'B': 18},
 # {'volume': 10, 't_s': 265, 't': 10, 'A': 13, 'B': 27},
 # {'volume': 10, 't_s': 270, 't': 10, 'A': 13, 'B': 31},
 # {'volume': 10, 't_s': 275, 't': 10, 'A': 13, 'B': 32},
 # {'volume': 10, 't_s': 280, 't': 10, 'A': 13, 'B': 33},
 # {'volume': 10, 't_s': 285, 't': 10, 'A': 13, 'B': 35},
 # {'volume': 10, 't_s': 290, 't': 10, 'A': 15, 'B': 36},
 # {'volume': 10, 't_s': 295, 't': 10, 'A': 16, 'B': 28},
 # {'volume': 10, 't_s': 300, 't': 10, 'A': 16, 'B': 34},
 # {'volume': 10, 't_s': 305, 't': 10, 'A': 17, 'B': 27}]
 # {'volume': 10, 't_s': 310, 't': 10, 'A': 17, 'B': 31},
 # {'volume': 10, 't_s': 315, 't': 10, 'A': 17, 'B': 32},
 # {'volume': 10, 't_s': 320, 't': 10, 'A': 17, 'B': 33},
 # {'volume': 10, 't_s': 325, 't': 10, 'A': 17, 'B': 35},
 # {'volume': 10, 't_s': 330, 't': 10, 'A': 18, 'B': 19},
 # {'volume': 10, 't_s': 335, 't': 10, 'A': 18, 'B': 20},
 # {'volume': 10, 't_s': 340, 't': 10, 'A': 18, 'B': 28},
 # {'volume': 10, 't_s': 345, 't': 10, 'A': 18, 'B': 34}]
 # {'volume': 10, 't_s': 350, 't': 10, 'A': 19, 'B': 36},
 # {'volume': 10, 't_s': 355, 't': 10, 'A': 20, 'B': 36},
 # {'volume': 10, 't_s': 360, 't': 10, 'A': 23, 'B': 27},
 # {'volume': 10, 't_s': 365, 't': 10, 'A': 23, 'B': 31},
 # {'volume': 10, 't_s': 370, 't': 10, 'A': 23, 'B': 32},
 # {'volume': 10, 't_s': 375, 't': 10, 'A': 23, 'B': 33},
 # {'volume': 10, 't_s': 380, 't': 10, 'A': 23, 'B': 35},
 # {'volume': 10, 't_s': 385, 't': 10, 'A': 24, 'B': 27}]
 # {'volume': 10, 't_s': 390, 't': 10, 'A': 24, 'B': 31},
 # {'volume': 10, 't_s': 395, 't': 10, 'A': 24, 'B': 32},
 # {'volume': 10, 't_s': 400, 't': 10, 'A': 24, 'B': 33},
 # {'volume': 10, 't_s': 405, 't': 10, 'A': 24, 'B': 35},
 # {'volume': 10, 't_s': 410, 't': 10, 'A': 25, 'B': 36},
 # {'volume': 10, 't_s': 415, 't': 10, 'A': 26, 'B': 27},
 # {'volume': 10, 't_s': 420, 't': 10, 'A': 26, 'B': 31},
 # {'volume': 10, 't_s': 425, 't': 10, 'A': 26, 'B': 32},
 # {'volume': 10, 't_s': 430, 't': 10, 'A': 26, 'B': 33},
 # {'volume': 10, 't_s': 435, 't': 10, 'A': 26, 'B': 35}]
 # {'volume': 10, 't_s': 440, 't': 10, 'A': 0, 'B': 16},
 # {'volume': 10, 't_s': 445, 't': 10, 'A': 0, 'B': 18},
 # {'volume': 10, 't_s': 450, 't': 10, 'A': 1, 'B': 16},
 # {'volume': 10, 't_s': 455, 't': 10, 'A': 1, 'B': 18},
 # {'volume': 10, 't_s': 460, 't': 10, 'A': 1, 'B': 36},
 # {'volume': 10, 't_s': 465, 't': 10, 'A': 2, 'B': 16},
 # {'volume': 10, 't_s': 470, 't': 10, 'A': 2, 'B': 18},
 # {'volume': 10, 't_s': 475, 't': 10, 'A': 4, 'B': 9},
 # {'volume': 10, 't_s': 480, 't': 10, 'A': 4, 'B': 10},
 # {'volume': 10, 't_s': 485, 't': 10, 'A': 4, 'B': 12},
 # {'volume': 10, 't_s': 490, 't': 10, 'A': 4, 'B': 13},
 # {'volume': 10, 't_s': 495, 't': 10, 'A': 4, 'B': 17},
 # {'volume': 10, 't_s': 500, 't': 10, 'A': 4, 'B': 23},
 # {'volume': 10, 't_s': 505, 't': 10, 'A': 4, 'B': 24},
 # {'volume': 10, 't_s': 510, 't': 10, 'A': 4, 'B': 26},
 # {'volume': 10, 't_s': 515, 't': 10, 'A': 4, 'B': 27},
 # {'volume': 10, 't_s': 520, 't': 10, 'A': 4, 'B': 31},
 # {'volume': 10, 't_s': 525, 't': 10, 'A': 4, 'B': 32},
 # {'volume': 10, 't_s': 530, 't': 10, 'A': 4, 'B': 33},
 # {'volume': 10, 't_s': 535, 't': 10, 'A': 4, 'B': 35},
 # {'volume': 10, 't_s': 540, 't': 10, 'A': 5, 'B': 9},
 # {'volume': 10, 't_s': 545, 't': 10, 'A': 5, 'B': 10},
 # {'volume': 10, 't_s': 550, 't': 10, 'A': 5, 'B': 12},
 # {'volume': 10, 't_s': 555, 't': 10, 'A': 5, 'B': 13},
 # {'volume': 10, 't_s': 560, 't': 10, 'A': 5, 'B': 17},
 # {'volume': 10, 't_s': 565, 't': 10, 'A': 5, 'B': 23},
 # {'volume': 10, 't_s': 570, 't': 10, 'A': 5, 'B': 24},
 # {'volume': 10, 't_s': 575, 't': 10, 'A': 5, 'B': 26},
 # {'volume': 10, 't_s': 580, 't': 10, 'A': 5, 'B': 27},
 # {'volume': 10, 't_s': 585, 't': 10, 'A': 5, 'B': 31},
 # {'volume': 10, 't_s': 590, 't': 10, 'A': 5, 'B': 32},
 # {'volume': 10, 't_s': 595, 't': 10, 'A': 5, 'B': 33}]

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
print("Timeline ")
print(timeline)
print("Time intervals ")
print(Intervals)

print("MIP evaluationg ...")
status, correct_mip, answ_mip, time_mip = mip.evaluate_model(V, K, M, Contracts, Intervals, Bandwidth, BackFlows)
print("Total time: ", time_mip)
print_answer(status, correct_mip, answ_mip, V, K, M)

print("Pyomo evaluationg ...")
status, correct_pyomo, answ_pyomo, time_pyomo = spyomo.evaluate_model(V, K, M, Contracts, Intervals, Bandwidth, BackFlows)
print("Total time: ", time_pyomo)
print_answer(status, correct_pyomo, answ_pyomo, V, K, M)

print("Does MIP and Pyomo answers are equal?")
print((answ_mip == answ_pyomo).all())


# In[ ]:




