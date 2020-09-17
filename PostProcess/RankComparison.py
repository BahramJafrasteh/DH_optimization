import scipy.io as spio
import os.path
import numpy as np
from pathlib import Path
from collections import defaultdict
from numpy import loadtxt
from copy import copy
from sys import platform
import pickle as pickle
from scipy.stats import norm #pdf and cdf
from random import randint # random integer
def main():
    ESF_data = True
    SAR_data = False
    if platform == 'win32':
        if ESF_data:
            Strs = 'Esfordi'
            path_read_predicted = 'C:/Model/'
            path_read = 'C:/CV'
            path_read_part = 'C:/RandBH_ang_877/'
        elif SAR_data:
            Strs = 'Sarcheshmeh'
            path_read_predicted = 'C:/Predict/'
            path_read = 'C:/Sarcheshmeh_CV/CV'
            path_read_part = 'C://RandBH_ang/'
    else:
        if ESF_data:
            Strs = 'Esfordi'
            path_read_predicted = '/Esfordi_GP_Symmetric_ExpoAns_final/Model/'
            path_read = '/CV'
            path_read_part = '/RandBH_ang_877/'
        elif SAR_data:
            Strs = 'Sarcheshmeh'
            path_read_predicted = '/Sarcheshmeh_CV/Model/'
            path_read = '/Sarcheshmeh_CV/CV'
            path_read_part = '/RandBH_ang/'


    rr = 0

    nums = 1
    maxrep = 50
    maxnumadd = 17#196#17
    #196

    rm = 0
    nrep = 0
    Kendals = []
    Footrules = []
    while nrep < maxrep:
        nrep = nrep + 1
        kl = 50
        Kendal_vals = np.zeros(shape = (maxnumadd,7), dtype = float)
        Footrule_vals = np.zeros(shape = (maxnumadd,7), dtype = float)
        for k in range(0, kl):
            counter = 0
            fname = path_read_part
            fname += 'data_divide_NN' + str(k+1) + '.mat'
            pathfile = Path(fname)
            if pathfile.is_file():
                #print("rading file ...")
                rdata = spio.loadmat(fname, squeeze_me=True)
                #print("done")
            else:
                print("Error file does not exist")
            if ESF_data:
                ALPHA_GET = 2.0
                BETA_GET = 0.22
                GAMMA_GET = 0.22
                ALPHA_GE = 2.0
                BETA_GE = 2.0
                ALPHA_BJ = 1.44
                BETA_BJ = 1.55
                GAMMA_BJ = 1.55
                thrsh_pi = 0.7
                thrsh_ei = 0.8
                fpredt = path_read_predicted + 'Esfordi_test_predict_' + str(k+1) + '_Desort.txt'

            elif SAR_data:
                ALPHA_GET = 1.88
                BETA_GET = 0.88
                GAMMA_GET = 0.88
                ALPHA_GE = 2.0
                BETA_GE = 2.0
                ALPHA_BJ = 2.0
                BETA_BJ = 0.22
                GAMMA_BJ = 0.22
                thrsh_pi = 0.1
                thrsh_ei = 0.2
                fpredt = path_read_predicted + 'Sarch_test_' + str(k+1) + '_ExpAns_predict_Desort.txt'

            pathfile = Path(fpredt)
            if pathfile.is_file():
                #print("rading file ...")
                Mt = loadtxt(fpredt, comments="#", delimiter="\t", unpack=False)
                #print("done")
            else:
                print("Error file does not exist")
            Ytr = rdata['Ytr']
            Yst = Mt[:,3]
            stdYst = max(Yst) - min(Yst)
            meanYst = min(Yst)
            Yht = Mt[:,2]
            Yht[ Yht < 0.0 ] = 0.0
            meanYh = min(Yht)
            stdYh = max(Yht) - min(Yht)
            if ESF_data:
                Rock = rdata['Xtr'][:,4]
                Rock_test = rdata['Xt'][:,4]
            elif SAR_data:
                Rock = rdata['Xtr'][:,5]
                Rock_test = rdata['Xt'][:,5]
            BHID_t = rdata['Xt'][:,0]
            U_BH = np.unique(BHID_t)
            Ytr_norm = ( Ytr-min(Ytr) )/ ( max(Ytr)- min(Ytr)  )

            L_rock = []
            U_rock = np.unique(Rock)
            for rn in range(0, len(U_rock)):
                ind_rck = (U_rock[rn] == Rock)
                L_rock.append( np.sum(Ytr_norm[ind_rck]) )
            L_rock = 1 - L_rock / np.sum(Ytr_norm)
            yobserved = np.zeros(shape=(len(U_BH), ), dtype = float)
            Func_GET = np.zeros(shape=(len(U_BH), ), dtype = float)
            Func_EEV = np.zeros(shape=(len(U_BH), ), dtype = float)
            Func_GE = np.zeros(shape=(len(U_BH), ), dtype = float)
            Func_BJ = np.zeros(shape=(len(U_BH), ), dtype = float)
            Func_PI = np.zeros(shape=(len(U_BH), ), dtype = float)
            Func_EI = np.zeros(shape=(len(U_BH), ), dtype = float)
            Thick = np.empty(shape=(len(U_BH), ), dtype = float)
            Yht_norm = ( Yht - min(Yht) )  / (max(Yht) - min(Yht) )
            Yst_norm = ( Yst - min(Yst) )  / (max(Yst) - min(Yst) )
            for i in range(0, len(U_BH)):
                ind_BH = BHID_t == U_BH[i]
                Yht_use = Yht[ind_BH];
                Thick[i] = len(Yht_use)
            Thick = (Thick - min(Thick))/ (max(Thick) - min(Thick))
            for i in range(0, len(U_BH)):
                #index of the ith drillhole
                ind_BH = BHID_t == U_BH[i];
                Yht_use = Yht[ind_BH]
                Yst_use = Yst[ind_BH]
                Rock_test_use = Rock_test[ind_BH]
                phosphate_price = 100 #$/ton
                density_phosphate = 2.39 #ton/m3
                V_block = 1000 #1000m3
                Processing_p = 1.1 #$/ton
                recovery_rate = 0.9 #*100
                MinigCost = 4.0 #$/ton
                cut_off = 5.0
                yobserved[i] = np.median( rdata['Yt'][ind_BH] )#*len( rdata['Yt'][ind_BH] )#np.mean( rdata['Yt'][ind_BH] )
                maxYh = max(Yht)
                if (maxYh != 0 and max(Yht) != min(Yht) ):
                    Yht_BO = ( Yht_use - min(Yht) )  / (maxYh - min(Yht) )
                else:
                    Yht_BO = ( Yht_use - min(Yht) )
                maxYst = max(Yst)
                if (maxYst != 0):
                    Yst_BO = ( Yst_use - min(Yst) ) /( maxYst - min(Yst) )
                else:
                    Yst_BO = ( Yst_use - min(Yst) )

                Func_GET[i] = np.dot(pow(Yht_BO,ALPHA_GET), pow(Yst_BO,BETA_GET))  * pow(Thick[i],GAMMA_GET)

                Func_GE[i] = np.dot( pow(Yht_BO, ALPHA_GE) , pow(Yst_BO, BETA_GE) )
                U_rock_test = np.unique(Rock_test_use)
                tt = np.empty(shape=(len(Yht_BO), ), dtype = float)
                ind_rck = (U_rock_test[0] == Rock_test_use)
                tt[ind_rck] = np.multiply( pow(Yht_BO[ind_rck]  ,ALPHA_BJ)
                                         , pow(Yst_BO[ind_rck], BETA_BJ) ) / ( pow(L_rock[0], GAMMA_BJ))
                rn = 1
                while (rn < len(U_rock_test)):
                    ind_rck = (U_rock_test[rn] == Rock_test_use);
                    tt[ind_rck] = np.multiply( pow(Yht_BO[ind_rck]  ,ALPHA_BJ)
                                         , pow(Yst_BO[ind_rck], BETA_BJ) ) / ( pow(L_rock[rn], GAMMA_BJ))
                    rn = rn + 1;
                Func_BJ[i] = np.mean(tt)
                # EEV
                nexy2 = Yht_use
                EEV = V_block*(density_phosphate *( phosphate_price - Processing_p) * recovery_rate * nexy2 - MinigCost )
                indth = np.argwhere( nexy2 < cut_off)
                EEV[indth] = -V_block*MinigCost* np.ones(shape=(len(indth),1), dtype=float)
                Func_EEV[i] =  np.mean(EEV)

                # Bayesian
                #Yht_BO = Yht_use
                #Yst_BO = Yst_use    #indl = Yht_BO < (cutoff - min(Yht)) / ( max(Yht) - min(Yht) );
                #thrsh_pi = thrsh_pi * max(Yht_BO)
                #thrsh_ei = thrsh_ei * max(Yht_BO)
                if ( np.median(Yht_BO) > 0 and np.median(Yst_BO) > 0 ):
                    Func_PI[i] = norm.cdf( np.median(Yht_BO),thrsh_pi ,np.median (Yst_BO) )
                    Func_EI[i]   = ( np.median(Yht_BO)-thrsh_ei ) * norm.cdf(np.median(Yht_BO),thrsh_ei,np.median(Yst_BO) ) +\
                    np.median(Yst_BO) * norm.pdf(np.median(Yht_BO),thrsh_ei,np.median(Yst_BO) )

                '''
                if (Yht_BO[0] > 0 and Yst_BO[0] > 0):
                    Func_PI[i] = norm.cdf(Yht_BO[0],thrsh_pi ,Yst_BO[0])
                    Func_EI[i]   = ( Yht_BO[0]-thrsh_ei ) * norm.cdf( Yht_BO[0],thrsh_ei, Yst_BO[0] ) +\
                        Yst_BO[0] * norm.pdf( Yht_BO[0],thrsh_ei, Yst_BO[0] )
                else:
                    Func_PI[i] = 0
                    Func_EI[i] = 0
                for ij in range(0, len(Yht_BO)):
                    if (Yht_BO[ij] > 0 and Yst_BO[ij] > 0):
                        Func_PI[i] = Func_PI[i] + norm.cdf(Yht_BO[ij],thrsh_pi ,Yst_BO[ij])
                        Func_EI[i]   = Func_PI[i] + ( Yht_BO[ij]-thrsh_ei ) * norm.cdf( Yht_BO[ij],thrsh_ei, Yst_BO[ij] ) +\
                            Yst_BO[ij] * norm.pdf( Yht_BO[ij],thrsh_ei, Yst_BO[ij] )
                        '''
                #Func_PI[i] = Func_PI[i] / len(Yht_BO)
                #Func_EI[i] = Func_EI[i] / len(Yht_BO)

            for ik in range(1, maxnumadd):

                rnad_gen = np.random.permutation(len(yobserved))
                rnad_gen = rnad_gen[range (0, ik+1)]
                #rnad_gen = np.array(list(range (0, ik+1)))
                random_val = sorted( np.random.permutation(ik+1), reverse = True)

                S_BH = sorted(yobserved[rnad_gen], reverse = True)
                S_BH_ind = np.array (list (reversed( np.argsort(yobserved[rnad_gen]) ) ) )
                S_GET = sorted(Func_GET[rnad_gen], reverse = True)
                S_GET_ind = np.array (list (reversed( np.argsort(Func_GET[rnad_gen]) ) ) )
                S_GE = sorted(Func_GE[rnad_gen], reverse = True)
                S_GE_ind = np.array (list (reversed( np.argsort(Func_GE[rnad_gen]) ) ) )
                S_BJ = sorted(Func_BJ[rnad_gen], reverse = True)
                S_BJ_ind = np.array (list (reversed( np.argsort(Func_BJ[rnad_gen]) ) ) )
                S_EEV = sorted(Func_EEV[rnad_gen], reverse = True)
                S_EEV_ind = np.array (list (reversed( np.argsort(Func_EEV[rnad_gen]) ) ) )
                S_PI = sorted(Func_PI[rnad_gen], reverse = True)
                S_PI_ind = np.array (list (reversed( np.argsort(Func_PI[rnad_gen]) ) ) )
                S_EI = sorted(Func_EI[rnad_gen], reverse = True)
                S_EI_ind = np.array (list (reversed( np.argsort(Func_EI[rnad_gen]) ) ) )
                if ( max(S_BH) == min (S_BH) ):
                    #weights = S_BH*1.0
                    weights = yobserved[rnad_gen]
                else:
                    #weights = (S_BH - min(S_BH) ) / ( max(S_BH)- min(S_BH) )
                    weights = (yobserved[rnad_gen] - min(S_BH) ) / ( max(S_BH)- min(S_BH) )
                tmp0 = copy(S_BH_ind)
                tmp1 = copy(S_GET_ind)
                tmp2 = copy(S_GE_ind)
                tmp3 = copy(S_BJ_ind)
                tmp4 = copy(S_PI_ind)
                tmp5 = copy(S_EI_ind)
                tmp6 = copy(S_EEV_ind)
                for ii in range(0, len(S_BH)):
                    indf = (S_GET_ind == S_BH_ind[ii])
                    tmp1[indf] = ii
                    indf = (S_GE_ind == S_BH_ind[ii])
                    tmp2[indf] = ii
                    indf = (S_BJ_ind == S_BH_ind[ii])
                    tmp3[indf] = ii
                    indf = (S_PI_ind == S_BH_ind[ii])
                    tmp4[indf] = ii
                    indf = (S_EI_ind == S_BH_ind[ii])
                    tmp5[indf] = ii
                    indf = (S_EEV_ind == S_BH_ind[ii])
                    tmp6[indf] = ii
                    tmp0[ii] = ii
                S_BH_ind = tmp0
                S_GET_ind = tmp1
                S_GE_ind = tmp2
                S_BJ_ind = tmp3
                S_PI_ind = tmp4
                S_EI_ind = tmp5
                S_EEV_ind = tmp6
                Kendal_vals[ik,0] = Kendal_vals[ik,0] + WeightedKendalTau(S_GET_ind, weights, ik)
                Kendal_vals[ik,1] = Kendal_vals[ik,1] + WeightedKendalTau(S_GE_ind, weights, ik)
                Kendal_vals[ik,2] = Kendal_vals[ik,2] + WeightedKendalTau(S_BJ_ind, weights, ik)
                Kendal_vals[ik,3] = Kendal_vals[ik,3] + WeightedKendalTau(S_PI_ind, weights, ik)
                Kendal_vals[ik,4] = Kendal_vals[ik,4] + WeightedKendalTau(S_EI_ind, weights, ik)
                Kendal_vals[ik,5] = Kendal_vals[ik,5] + WeightedKendalTau(random_val, weights, ik)
                Kendal_vals[ik,6] = Kendal_vals[ik,6] + WeightedKendalTau(S_EEV_ind, weights, ik)

                Footrule_vals[ik,0] = Footrule_vals[ik,0] + WeightedFootrule(S_GET_ind, weights, ik)
                Footrule_vals[ik,1] = Footrule_vals[ik,1] + WeightedFootrule(S_GE_ind, weights, ik)
                Footrule_vals[ik,2] = Footrule_vals[ik,2] + WeightedFootrule(S_BJ_ind, weights, ik)
                Footrule_vals[ik,3] = Footrule_vals[ik,3] + WeightedFootrule(S_PI_ind, weights, ik)
                Footrule_vals[ik,4] = Footrule_vals[ik,4] + WeightedFootrule(S_EI_ind, weights, ik)
                Footrule_vals[ik,5] = Footrule_vals[ik,5] + WeightedFootrule(random_val, weights, ik)
                Footrule_vals[ik,6] = Footrule_vals[ik,6] + WeightedFootrule(S_EEV_ind, weights, ik)
            print (Strs + ', repNo. :', str(nrep), ' numpar: ' + str(k) +
                        ', Wkendal_GET : ' +  str(Kendal_vals[ik,0]) +
                        ', Wkendal_GE : ' +  str(Kendal_vals[ik,1]) +
                        ', Wkendal_BJ : ' + str(Kendal_vals[ik,2]) +
                        ', Wkendal_PI : ' + str(Kendal_vals[ik,3]) +
                        ', Wkendal_EI : ' + str(Kendal_vals[ik,4]) +
                        ', Wkendal_Random : ' + str(Kendal_vals[ik,5])+
                        ', EEV : ' + str(Kendal_vals[ik,6])
                        )

        Footrule_vals = Footrule_vals / kl
        Kendal_vals = Kendal_vals / kl
        Kendals.append(copy(Kendal_vals))
        Footrules.append(copy(Footrule_vals))
    Kendal_vals.fill(0.0)
    Footrule_vals.fill(0.0)
    Std_Kendal = copy(Kendal_vals)
    Std_Footrule = copy(Footrule_vals)
    for i in range(0,len(Kendals)):
        Footrule_vals = Footrule_vals + np.array(Footrules[i])
        Kendal_vals = Kendal_vals + np.array(Kendals[i])
    Footrule_vals = Footrule_vals/len(Kendals)
    Kendal_vals = Kendal_vals/len(Kendals)

    for i in range(0,len(Kendals)):
        Std_Footrule = Std_Footrule + pow(np.array(Footrules[i]) - Footrule_vals, 2.0)
        Std_Kendal = Std_Kendal + pow( np.array(Kendals[i]) - Kendal_vals, 2.0)
    Std_Footrule = np.sqrt(Std_Footrule/ len(Kendals))
    Std_Kendal = np.sqrt(Std_Kendal/ len(Kendals))
    csv_columns = 'GET \t GE \t BJ \t PI \t EI \t Random\t EEV'
# Writing out results
    np.savetxt(Strs + '_Footrule.csv', Footrule_vals, delimiter = '\t', header = csv_columns)
    np.savetxt(Strs + '_Footrule_std.csv', Std_Footrule, delimiter = '\t', header = csv_columns)
    np.savetxt(Strs + 'Kendal.csv', Kendal_vals, delimiter = '\t', header = csv_columns)
    np.savetxt(Strs + 'Kendal_std.csv', Std_Kendal, delimiter = '\t', header = csv_columns)

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
