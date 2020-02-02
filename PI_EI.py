###
# Compute PI and EI values
###

import scipy.io as spio
import os.path
import numpy as np
from pathlib import Path
from collections import defaultdict
from numpy import loadtxt
from copy import copy
import pickle as pickle
from scipy.stats import norm #pdf and cdf
from random import randint # random integer
def main():
    ESF_data = False
    SAR_data = True
    if ESF_data:
        Strs = 'Esfordi'
        #path _predicted data
        path_read_predicted = ''

    elif SAR_data:
        Strs = 'Sarcheshmeh'
        # path predicted data
        path_read_predicted = ''


    rr = 0

    nums = 1
    sq1 = frange(0.0, 7.0, 0.05)
    maxnumadd = 15
    Taus_var = defaultdict(list)
    Taus_var['gamma'] = []
    Taus_var['NumberOfAddBHs'] = []
    Taus_var['PI_footr'] = []
    Taus_var['PI_kendal'] = []
    Taus_var['EI_footr'] = []
    Taus_var['EI_kendal'] = []

    rm = 0
    al = 0


    cl = 0
    while cl < len(sq1):
        ga =sq1[cl]
        cl = cl + 1
        thrsh_pi = ga
        thrsh_ei = ga
        for numadd in range(5, maxnumadd):
            kl = 50
            Kendal_vals = np.empty(shape = (kl,2), dtype = float)
            Footrule_vals = np.empty(shape = (kl,2), dtype = float)
            for k in range(0, kl):
                numcv = 2
                counter = 0
                WeightedFootruleDistance = np.empty(shape = (numcv,2), dtype = float)
                WKendalTauDistance = np.empty(shape = (numcv,2), dtype = float)
                for cv in range(0, numcv):
                    fname = path_read
                    if ESF_data:
                        fname += '/Esfordi_' + str(k+1) + '_cv_' + str(cv) + '.mat'
                    if SAR_data:
                        fname += '/Sarch_' + str(k+1) + '_cv_' + str(cv) + '.mat'
                    pathfile = Path(fname)
                    if pathfile.is_file():
                        #print("rading file ...")
                        rdata = spio.loadmat(fname, squeeze_me=True)
                        #print("done")
                    else:
                        print("Error file does not exist")
                    Ytr = rdata['Ytr']

                    if (ESF_data):
                        fpredt = path_read_predicted + 'Esfordi_test_' + str(k+1) +'_predict_cv_' + str(cv) + '_Desort.txt'
                    #Sarcheshmeh
                    elif(SAR_data):
                        fpredt = path_read_predicted +'Sarch_test_' + str(k+1) +'_predict_cv_'+ str(cv) + '_Desort.txt'
                    pathfile = Path(fpredt)
                    if pathfile.is_file():
                        #print("rading file ...")
                        Mt = loadtxt(fpredt, comments="#", delimiter="\t", unpack=False)
                        #print("done")
                    else:
                        print("Error file does not exist")
                    Yst = Mt[:,3]
                    Yht = Mt[:,2]
                    Yht[ Yht < 0.0 ] = 0.0
                    BHID_t = rdata['Xt'][:,0]
                    U_BH = np.unique(BHID_t)
                    yobserved = np.empty(shape=(len(U_BH), ), dtype = float)
                    Func_PI = np.empty(shape=(len(U_BH), ), dtype = float)
                    Func_EI = np.empty(shape=(len(U_BH), ), dtype = float)

                    for i in range(0, maxnumadd):
                        #index of the ith drillhole
                        ind_BH = BHID_t == U_BH[i];
                        Yht_use = Yht[ind_BH]
                        Yst_use = Yst[ind_BH]
                        yobserved[i] = np.mean( rdata['Yt'][ind_BH] )

                        Yht_BO = Yht_use
                        Yst_BO = Yst_use    #indl = Yht_BO < (cutoff - min(Yht)) / ( max(Yht) - min(Yht) );

                        # Bayesian
                        Func_PI[i] = norm.cdf( thrsh_pi ,np.mean(Yht_BO),np.mean (Yst_BO) )
                        Func_EI[i]   = ( thrsh_ei-np.mean(Yht_BO) ) * norm.cdf(thrsh_ei,np.mean(Yht_BO),np.mean(Yst_BO) ) +\
                            np.mean(Yst_BO) * norm.pdf(thrsh_ei,np.mean(Yht_BO),np.mean(Yst_BO) )
                    ik = numadd
                    S_BH = sorted(yobserved[range(0, ik)], reverse = True)
                    S_BH_ind = np.array (list (reversed( np.argsort(yobserved[range(0, ik)]) ) ) )
                    S_PI = sorted(Func_PI[range(0, ik)])
                    S_PI_ind = np.array (list (reversed( np.argsort(Func_PI[range(0, ik)]) ) ) )
                    S_EI = sorted(Func_EI[range(0, ik)])
                    S_EI_ind = np.array (list (reversed( np.argsort(Func_EI[range(0, ik)]) ) ) )
                    if ( max(S_BH) == min (S_BH) ):
                        weights = [1] * len(S_BH)
                    else:
                        weights = (S_BH - min(S_BH) ) / ( max(S_BH)- min(S_BH) )
                    tmp0 = copy(S_BH_ind)
                    tmp1 = copy(S_PI_ind)
                    tmp2 = copy(S_EI_ind)
                    for ii in range(0, len(S_BH)):
                        indf = (S_PI_ind == S_BH_ind[ii])
                        tmp1[indf] = ii
                        indf = (S_EI_ind == S_BH_ind[ii])
                        tmp2[indf] = ii
                        tmp0[ii] = ii
                    S_BH_ind = tmp0
                    S_PI_ind = tmp1
                    S_EI_ind = tmp2
                    WKendalTauDistance[cv,0] = WeightedKendalTau(S_PI_ind, weights, ik)
                    WKendalTauDistance[cv,1] = WeightedKendalTau(S_EI_ind, weights, ik)
                    WeightedFootruleDistance[cv,0] = WeightedFootrule(S_PI_ind, weights, ik)
                    WeightedFootruleDistance[cv,1] = WeightedFootrule(S_EI_ind, weights, ik)
                    print (Strs +  ', numadd : ' +  str(numadd) + ', numcv: ' + str(cv) +
                           ', numpar: ' + str(k) +
                            ', gamma: ' +  str(ga) +  ', Wkendal_PI : ' +  str(WKendalTauDistance[cv,0]) +  ', Wkendal_EI : ' +
                            str(WKendalTauDistance[cv,1])  )
                Kendal_vals[k, 0] = np.mean(WKendalTauDistance[:,0])
                Kendal_vals[k, 1] = np.mean(WKendalTauDistance[:,1])
                Footrule_vals[k, 0] = np.mean(WeightedFootruleDistance[:,0])
                Footrule_vals[k, 1] = np.mean(WeightedFootruleDistance[:,1])
            Taus_var['gamma'].append(ga)
            Taus_var['NumberOfAddBHs'].append(numadd)
            Taus_var['PI_footr'].append(np.mean(Footrule_vals[:,0]))
            Taus_var['PI_kendal'].append(np.mean(Kendal_vals[:, 0]))
            Taus_var['EI_footr'].append(np.mean(Footrule_vals[:,1]))
            Taus_var['EI_kendal'].append(np.mean(Kendal_vals[:,1]))

    # Writing out results
    fout = Strs + '_' + str(ga) + '_all.csv'
    data = np.empty(shape = (len(Taus_var['gamma']), len(Taus_var)), dtype = float)
    i = 0
    csv_columns = ''
    for keys in Taus_var:
        data[:,i] = Taus_var[keys]
        i = i + 1
        csv_columns += keys + ' \t '
    np.savetxt(fout, data, delimiter = '\t', header = csv_columns)

    Taus = defaultdict(list)
    Taus['gamma'] = []
    Taus['PI_F'] = []
    Taus['PI_K'] = []
    Taus['EI_F'] = []
    Taus['EI_K'] = []
    cl = 0
    while cl < len(sq1):
        ga =sq1[cl]
        cl = cl + 1
        ind = (data[:,2] == ga)
        dr = data[ind,:]
        Taus['gamma'].append(ga)
        Taus['PI_F'].append(np.mean(dr[:,2]))
        Taus['PI_K'].append(np.mean(dr[:,3]))
        Taus['EI_F'].append(np.mean(dr[:,4]))
        Taus['EI_K'].append(np.mean(dr[:,5]))

    indexPIK = np.argmin(Taus['PI_K'])
    indexEIK = np.argmin(Taus['PI_K'])
    indexPIF = np.argmin(Taus['PI_F'])
    indexEIF = np.argmin(Taus['PI_F'])
    d1 = [ Taus['gamma'][indexPIK], Taus['gamma'][indexPIF]  ]
    d2 = [Taus['gamma'][indexEIK], Taus['gamma'][indexEIF] ]
    fout = Strs + '_' + str(ga) + '_mean.csv'
    data2 = np.empty(shape = (len(Taus['gamma']), len(Taus)), dtype = float)
    i = 0
    csv_columns2 = ''
    for keys in Taus:
        data2[:,i] = Taus[keys]
        i = i + 1
        csv_columns2 += keys + ' \t '
    np.savetxt(fout, data2, delimiter = '\t', header = csv_columns2)

    print (['Objs','Kendal', 'Footrule'])
    print (['Optimal GET',d1])
    print (['Optimal GE',d2])


def frange(start, stop, step):
    tm = start
    bm = [start]
    while tm < stop:
        bm.append(tm+step)
        tm = tm + step
    return bm
def WeightedKendalTau(objval, weight, ik):
    tau = 0
    for rk in range(0, ik):
        lr = 0
        if rk > 0:
            while True:
                if objval[rk] < objval[lr]:
                    tau = tau + weight[rk]*weight[lr]
                lr += 1
                if lr > (rk - 1):
                    break
    return tau
def WeightedFootrule(objval, weight, ik):
    ds = 0
    for ij in range(0, ik):
        ds = ds + weight[ij] *abs ( sum(weight[range(0, ij+1 )]) -  sum(weight[range(0, objval[ij]+1)]) )
    return ds

def save_obj(obj, name ):
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open( name + '.pkl', 'rb') as f:
        return pickle.load(f)

if __name__ == "__main__": main() # it is important to have it to write functions in any order that you want
