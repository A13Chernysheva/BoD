#!/usr/bin/env python
# coding: utf-8

# In[1]:


from mip import *
from itertools import product
from parse import *
import time
import numpy as np


# ,,,,,,,,,,,,,, Обнуление тех значений, между которыми нет каналов ,,,,,,,,,,,,,,

def zero_nconnect(model, q, K, V, Bandwidth):
    for (k, i, j) in product(K, V, V):
        if Bandwidth.at[i, j] == 0:
            model += q[i][j][k] == 0


# ,,,,,,,,,,,,,, Условие 1 ,,,,,,,,,,,,,,

def eval_cnstr1(model, q, M, V, I, Bandwidth, BackFlows):
    for tm in M:
        for i in V:
            for j in V:
                if Bandwidth.at[i, j] != 0:
                    if I[tm]['Ctm'] != set():
                        model += xsum(q[i][j][k] for k in I[tm]['Ctm']) <= (int)(Bandwidth.at[i, j] - BackFlows.at[i, j])


# ,,,,,,,,,,,,,, Условие 2 ,,,,,,,,,,,,,,

def eval_cnstr2(model, q, C, K, V, M, I):
    for k in K:
        a_k = C[k]['A']
        s = 0
        for tm in M:
            if k in I[tm]['Ctm']:
                for j in V:
                    s += q[a_k][j][k] * I[tm]['em']
        model += (s == C[k]['volume'])


# ,,,,,,,,,,,,,, Условие 3 ,,,,,,,,,,,,,,

def eval_cnstr3(model, q, C, K, V, M, I):
    for k in K:
        a_k = C[k]['A']
        s = 0
        for tm in M:
            if k in I[tm]['Ctm']:
                for j in V:
                    s += q[j][a_k][k] * I[tm]['em']
        model += (s == 0)


# ,,,,,,,,,,,,,, Условие 4 ,,,,,,,,,,,,,,

def eval_cnstr4(model, q, C, K, V, M, I):
    for k in K:
        b_k = C[k]['B']
        s = 0
        for tm in M:
            if k in I[tm]['Ctm']:
                for j in V:
                    s += q[j][b_k][k] * I[tm]['em']
        model += (s == C[k]['volume'])


# ,,,,,,,,,,,,,, Условие 5 ,,,,,,,,,,,,,,

def eval_cnstr5(model, q, C, K, V, M, I):
    for k in K:
        b_k = C[k]['B']
        s = 0
        for tm in M:
            if k in I[tm]['Ctm']:
                for j in V:
                    s += q[b_k][j][k] * I[tm]['em']
        model += (s == 0)


# ,,,,,,,,,,,,,, Условие 6 ,,,,,,,,,,,,,,

def eval_cnstr6(model, q, M, K, C, V, Bandwidth, I):
    for tm in M:
        for k in I[tm]['Ctm']:
            a_k = C[k]['A']
            b_k = C[k]['B']
            for r in V - {a_k, b_k}:
                sum_in = 0
                sum_out = 0
                flag = False
                for j1 in V:
                    if Bandwidth.at[r, j1] != 0:
                        sum_out += q[r][j1][k]
                        flag = True
                for j2 in V:
                    if Bandwidth.at[j2, r] != 0:
                        sum_in += q[j2][r][k]
                        flag = True
                if flag:
                    model += ((sum_in - sum_out) == 0)


# ,,,,,,,,,,,,,, Условие 7 ,,,,,,,,,,,,,,

def eval_cnstr7(model, q, V, K):
    for i in V:
        for j in V:
            for k in K:
                model += q[i][j][k] >= 0


def get_result(status, model, V, K, M, time):
    correctness = ""
    answer = np.zeros([len(V), len(V), len(K)])
    if (status != None) and (status == OptimizationStatus.OPTIMAL):
        count = 0
        for i in range(len(model.vars)):
            pattern = 'q[{:d}][{:d}][{:d}]'
            i1, i2, i3 = parse(pattern, model.vars[i].name)
            answer[i1][i2][i3] = model.vars[i].x
            if (answer[i1][i2][i3] != 0):
                count += 1 
        print('Number of non zero variables: ' + str(count))
        print('----------------------')
    else:
        correctness = "Model infeasible"
        
    return (len(model.constrs), correctness, answer, time)


def evaluate_model(V, K, M, C, I, Bandwidth, BackFlows):
    
    start_time = time.time()

    model = Model()

    q = [[[model.add_var(name = 'q[{}][{}][{}]'.format(i,j,k),
              var_type = INTEGER)  for k in K]
                                   for j in V] 
                                   for i in V]

    zero_nconnect(model, q, K, V, Bandwidth)
    eval_cnstr1(model, q, M, V, I, Bandwidth, BackFlows)
    eval_cnstr2(model, q, C, K, V, M, I)
    eval_cnstr3(model, q, C, K, V, M, I)
    eval_cnstr4(model, q, C, K, V, M, I)
    eval_cnstr5(model, q, C, K, V, M, I)
    eval_cnstr6(model, q, M, K, C, V, Bandwidth, I)
    eval_cnstr7(model, q, V, K)
    
    model.objective = minimize(xsum(q[i][j][k] for i, j, k in product(V, V, K)))

    print('----------------------')
    print('Total number of variables: ' + str(len(model.vars)))
    
    status = model.optimize()
    
    end_time = time.time()
    
    return get_result(status, model, V, K, M, end_time - start_time)



