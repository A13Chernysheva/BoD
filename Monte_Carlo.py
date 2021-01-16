#!/usr/bin/env python
# coding: utf-8

# Фиксируем все параметры, кроме времени наступления контрактов и расчитаем в скольки случаях подобный контракт будет удовлетворен.

# Для этого будем читать из файла *settings_mc.ini*:
# * имя топологии
# * характеристики контрактов, чтобы найти соответствующий файл в Unique_Contracts
# * количество проводимых исптыаний
# 

# Основные шаги:
# 1. Считать топологию по имени, в ней уже зафиксированы фоновые потоки и пропускные способности
# 2. Считать уникальный набор контрактов
# 3. Проводим N шагов метода Монте-Карло, где на каждом шаге:\
#     a. Формуруем финальный набор контрактов (фиксируем время старта)\
#     b. Вычисляем результат\
#     c. Запоминаем смогли бы мы удовлетворить данный контракт или нет
# 4. Показываем результат - массив из 0 и 1, где позиция числа - номер итерации, 0 - не смогли, 1 - смогли

# In[16]:


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

from tqdm.notebook import tqdm

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

def def_start_time(Unique_contracts, ACON_NUM, UCON_NUM, T_MAX):
    Contracts = []
    for i in range(ACON_NUM):
        t_s = int(random.uniform(0, T_MAX))
        while(True):
            i = random.randint(0, UCON_NUM - 1)
            c = Unique_contracts[i]
            if type(c['t_s']) == int:
                d = [c['t_s']]
            else:
                d = c['t_s']
            if T_MAX - c['t'] < t_s:
                continue
            status = True
            for j in d:
                if t_s >= j and t_s <= j + c['t']:
                    status = False
            if not status:
                continue
            Unique_contracts[i]['t_s'] += [t_s]
            c = {'volume': Unique_contracts[i]['volume'],
                    't_s': t_s,
                      't': Unique_contracts[i]['t'],
                      'A': Unique_contracts[i]['A'],
                      'B': Unique_contracts[i]['B'],
                      'type': Unique_contracts[i]['type']}
            Contracts += [c]
            break  
    return Contracts

random.seed(3)

config = configparser.ConfigParser()
config.read("settings_mc.ini")
TOPO = config["Def"]["TOPO_NAME"]
CTYPE1 = int(config["Def"]["CTYPE1"])
CTYPE2 = int(config["Def"]["CTYPE2"])
CTYPE3 = int(config["Def"]["CTYPE3"])
CONTRACTS_NUM = int(config["Def"]["ACON_NUM"])
UCONTRACTS_NUM = CTYPE1 + CTYPE2 + CTYPE3
T_MAX = int(config["Def"]["T_MAX"])
N = int(config["Def"]["N"])

#  ----- Подготовка топологии -----

H1=nx.read_graphml('topo/{}.graphml'.format(TOPO))
H = nx.Graph(H1)

print("\nInfo about topology\n")
print("Name:   " + "{}".format(TOPO))
print("Number of nodes: " + str(len(H.nodes())))
print("Number of edges: " + str(len(H.edges())))

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
    
# ----- Считывание уникальных контрактов -----

print('\nInfo about contracts\n')
print('Number of contracts at all:', CONTRACTS_NUM)
print('Number of unique contracts:', UCONTRACTS_NUM)
print('Number of contracts TYPE 1:', CTYPE1)
print('Number of contracts TYPE 2:', CTYPE2)
print('Number of contracts TYPE 3:', CTYPE3)
print('T max:', T_MAX)

PATH = 'Unique_Contracts/{}'.format(TOPO)

with open('{}/C={}_CT1-{}_CT2-{}_CT3-{}.yaml'.format(PATH, UCONTRACTS_NUM, CTYPE1, CTYPE2, CTYPE3), 'r') as f:
	UContracts = yaml.load(f, Loader=yaml.FullLoader)

# ----- Начало метода Монте-Карло -----
    
print('********************************')
print("Number of iterations Monte-Carlo:", N)

Result = []

for n in tqdm(range(N)):
    Contracts = def_start_time(UContracts, CONTRACTS_NUM, UCONTRACTS_NUM, T_MAX)

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

    # print("Timeline len " + str(len(timeline)))
    # print(timeline)
    # print("Time intervals len " + str(len(Intervals)))
    # print(Intervals)

    print("MIP evaluating ...")
    status, correct_mip, answ_mip, time_mip = mip.evaluate_model(V, K, M, Contracts, Intervals, Bandwidth, BackFlows)
    #print(correct_mip)
    #print_answer(status, correct_mip, answ_mip, V, K)
    #print("Total time: ", time_mip)
    if correct_mip == 'correct':
        print("Pass")
        Result += ['1']
    else:
        print("Fail")
        Result += ['0']
    
print(Result)


# In[ ]:




