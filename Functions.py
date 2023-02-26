import numpy as np
import pandas as pd
import ROOT 
from math import sqrt


def Create_histo(column, name, bins, x0,x1):
    h = ROOT.TH1F(name, "",bins,x0,x1)
    for i in column:
        h.Fill(i)
    return h

def Create_histo_variable(I_array, name, bins, threshold =5, thr_problem=50, verbose=0):
    I_rms= I_array.std()
    I_mean= I_array.mean()
    I_min = I_mean-I_rms
    I_max = I_mean+I_rms
    h = ROOT.TH1F(name, "",bins,I_min,I_max)
    n_prob=0
    I_array_new = []
    if(verbose>0):
        print("h", name, " I_mean ", I_mean ," rms ", I_rms, " min ",I_min," max ",I_max)
    for i in I_array:
        if (abs(i-I_mean))>thr_problem:
            if (verbose>1):
                print("warning! Anomalous event with i= ",i," we are excluding it for now, it is the problematic event # ",n_prob," in this run")
            n_prob= n_prob+1
            continue
        if (abs(i-I_mean)/I_rms)<threshold:
            I_array_new.append(i) 
    I_array_new = np.array(I_array_new)
    I_Mean = I_array_new.mean()
    I_Rms = I_array_new.std()
    I_Min = I_Mean-I_Rms
    I_Max = I_Mean+I_Rms
    h = Create_histo(column = I_array_new,name = name, bins = 45, x0=I_Min, x1 =I_Max)

    return h, I_array_new
            
def Do_histo(folder, fname, f1,colum , col_del = "I3B",lim = 100, canvas= ROOT.TCanvas(), num = 1):
    parI = []
    sigI = []
    par_err = []
    sig_err = []
    num = len(fname)
    h = list(np.zeros(num))
    
    for i in range(num):
        df_I = pd.read_csv(folder + "/" + fname[i])
        if lim >0:
            df_I = df_I[df_I[col_del] <= lim]
        else:
            df_I = df_I[df_I[col_del] >= lim]
        I_arr = np.array(df_I[colum])
        canvas.cd(i+1)
        #if num == 1:
        h[i],I_arr_new= Create_histo_variable(I_arr, colum + " run " + str(i+1), 45)
        f1.SetParameters(1000,I_arr_new.mean(),I_arr_new.std())
        #else:
         #   h[i] = Create_histo(I_arr, colum + " run " + str(i+1), 45, I_arr.min(), I_arr.max())
          #  f1.SetParameters(1000,I_arr.mean(),I_arr.std())
        h[i].Fit(f1, "Q")
        parI.append(f1.GetParameter(1))
        sigI.append(abs(f1.GetParameter(2)))
        par_err.append(f1.GetParError(1))
        sig_err.append(f1.GetParError(2))
        

    return h,parI,sigI, par_err,sig_err

def V_mean(folder, fname, column):
    parV = []
    sigV = []
    num = len(fname)
    for i in range(num):
        df_V = pd.read_csv(folder + "/" +fname[i])
        col1 = df_V.columns[column]
        V1 = np.array(df_V[col1])
        if column != 0:
            col0 = df_V.columns[column-1]
            V0  = np.array(df_V[col0])
        else:
            V0 = np.zeros(len(V1))
        dV = V1-V0
        nV=np.trim_zeros(np.array(dV))
        meanV=nV.mean()
        stdV=nV.std()
        parV.append(meanV)
        sigV.append(abs(stdV))
    return parV,sigV

def K_mean(folder,fname):
    parIx = []
    parIy = []
    sigIx = []
    sigIy = []
    num = len(fname)
    for i in range(num):
        df_K = pd.read_csv(folder + "/" + fname[i])
        Ix = np.array(df_K["Ix"])
        Iy = np.array(df_K["Iy"])
        meanIx = Ix.mean()
        meanIy = Iy.mean()
        stdIx = Ix.std()
        stdIy = Iy.std()
        parIx.append(meanIx)
        parIy.append(meanIy)
        sigIx.append(abs(stdIx))
        sigIy.append(abs(stdIy))
    return parIx, parIy, sigIx, sigIy

def Do_histo_sum(folder, fname, f1,colum , col_del = "I3B",lim = 100, canvas= ROOT.TCanvas()):
    parI = []
    sigI = []
    par_err = []
    sig_err = []
    num = len(fname)
    h = list(np.zeros(num))
    
    for i in range(num):
        df_I = pd.read_csv(folder + "/" + fname[i])
        if lim >0:
            df_I = df_I[df_I[col_del] <= lim]
        else:
            df_I = df_I[df_I[col_del] >= lim]
        col = df_I.columns[colum]
        col0 = df_I.columns[colum -1]
        I_arr = np.array(df_I[col]) + np.array(df_I[col0])
        canvas.cd(i+1)
#       h[i] = Create_histo(I_arr, colum + " run " + str(i+1), 45, I_arr.min(), I_arr.max())
        h[i],I_arr_new= Create_histo_variable(I_arr, "sum run " + str(i+1), 45)
        f1.SetParameters(1000,I_arr_new.mean(),I_arr_new.std())
        h[i].Fit(f1, "Q")
        parI.append(f1.GetParameter(1))
        sigI.append(f1.GetParameter(2))
        par_err.append(f1.GetParError(1))
        sig_err.append(f1.GetParError(2))

    return h,parI,sigI, par_err, sig_err

def Do_histo_sum_total(folder, fname,f1,start, stop, col_del = "I3B",lim = 100, canvas= ROOT.TCanvas()):
    parI = []
    sigI = []
    par_err = []
    sig_err = []
    num = len(fname)
    h = list(np.zeros(num))
    
    for i in range(num):
        df_I = pd.read_csv(folder + "/" + fname[i])
        #df_K = pd.read_csv("arrays/kpa_data/" + fname2[i])
        if lim >0:
            df_I = df_I[df_I[col_del] <= lim]
        else:
            df_I = df_I[df_I[col_del] >= lim]
        col = df_I.columns[start]
        I_arr = np.array(df_I[col]) 
        for j in range(start+1,stop+1):
            col= df_I.columns[j]
            I_arr += np.array(df_I[col])
        canvas.cd(i+1)
        h[i] = Create_histo(I_arr, " run " + str(i+1), 45, I_arr.min(), I_arr.max())
        #h[i],I_arr_new= Create_histo_variable(I_arr, "sum run " + str(i+1), 45)
        f1.SetParameters(1000,I_arr.mean(),I_arr.std())
        h[i].Fit(f1, "Q")
        parI.append(f1.GetParameter(1))
        sigI.append(f1.GetParameter(2))
        par_err.append(f1.GetParError(1))
        sig_err.append(f1.GetParError(2))

    return h,parI,sigI, par_err, sig_err

