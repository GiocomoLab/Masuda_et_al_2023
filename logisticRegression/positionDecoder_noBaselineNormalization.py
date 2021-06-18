#!/usr/bin/env python
# coding: utf-8

# position_decoder.ipynb
# 
# Script to create a position decoder using multinomial logistic regression
# 
# MGC 7/15/2019
# FKM 7/31/2019

# In[171]:


# packages
import sys,os
#import pandas
import math
import glob
import pickle
from scipy import io
import numpy as np
from numpy import squeeze
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import loadmat as lm
from scipy.interpolate import make_interp_spline
from scipy.signal import savgol_filter
import time

# data
# neuropix_folder = os.path.join('E:\\','Dropbox','Work','neuropixels')
# data_dir = os.path.join(neuropix_folder,'data')
# sessions = pandas.read_excel(io=os.path.join(neuropix_folder,'sessions_range_of_contrasts.xlsx'))

# G4_190620_keicontrasttrack_baseline+cntrlinjx+ketamine.mat

def trainTestLogisticRegression(filename, params):

    # params
    
    dt = 0.2 # time bin for decoding, in ms
    every_nth_time_bin = round(dt/float(params['TimeBin']));
    dx = 5 # spatial bin for decoding, in cm
    track_start = float(params['TrackStart'])
    track_end = float(params['TrackEnd'])
    track_length = track_end-track_start
    numposbins = math.floor((track_end-track_start)/dx)
    posx_edges = np.linspace(track_start,track_end,numposbins+1)


    dat = lm.loadmat(filename)
    post = dat['post']
    posx = dat['posx']
    trial = dat['trial']
    sp = dat['sp']
    good_cells = sp['cids'][sp['cgs']==2]
    # spike_depth = dat['spike_depth']

    # sort cells by depth on probe
    # good_cells = good_cells[np.argsort(spike_depth)][::-1]
    # spike_depth = np.sort(spike_depth)[::-1]

    # resample post, posx, and trial according to desired dt
    post = post[0::every_nth_time_bin]
    posx = posx[0::every_nth_time_bin]
    trial = trial[0::every_nth_time_bin]

    # time bins for position decoding
    numtimebins = len(post)
    post_edges = squeeze(np.linspace(min(post)-dt/2,max(post)+dt/2,numtimebins+1))
    post_centers = post_edges[range(0,len(post_edges)-1)]+dt/2

    # posx categories for position decoding (binned)
    posx_bin = np.digitize(posx,posx_edges)

    # count spikes in each time bin for each cell
    spikecount = np.empty((len(good_cells),len(post),))
    spikecount[:] = np.nan
    for cell_idx in range(len(good_cells)):   
        spike_t = sp['st'][sp['clu']==good_cells[cell_idx]]
        spikecount[cell_idx,:] = np.histogram(spike_t,bins=post_edges)[0]


    # train the logistic regression model
    baseline_rmse = []
    train_set = np.squeeze((trial>=1) & (trial<=50))
    X_train = np.transpose(spikecount[:,train_set])
    y_train = squeeze(posx_bin[train_set])

    C_param = 0.03

    # train the logistic regression model 50 times dropping 1 to find the baseline RMSE
    # baseline_rmse = []
    # for x in range(1,51):
    #     try:
    #         train_set = np.squeeze((trial>=1) & (trial<=50) & (trial!=x))
    #         X_train = np.transpose(spikecount[:,train_set])
    #         y_train = squeeze(posx_bin[train_set])
    #         model = linear_model.LogisticRegression(penalty='l2',random_state=0, solver='lbfgs', multi_class='multinomial', max_iter=10000000, C = C_param).fit(X_train, y_train) 
    #         rmse = testModelRMSE(x,trial,spikecount,posx_bin,model,track_length,dx)
    #         baseline_rmse.append(rmse)
    #         sys.stdout.write('\n'+filename+'_baseline training: '+str(x))
    #         sys.stdout.flush()
    #     except Exception as e:
    #         sys.stdout.write('Could not train on: '+str(x))
    #         sys.stdout.write(str(e))
    #         sys.stdout.flush()


    # avgBaselineRMSE = np.mean(baseline_rmse)
    

    model = linear_model.LogisticRegression(penalty='l2',random_state=0, solver='lbfgs', multi_class='multinomial', max_iter=10000000, C = C_param).fit(X_train, y_train)

    all_rmse = []
    plotRange = range(50, 300, 1)
    for x in plotRange:
        # sys.stdout.write(filename,': ',x)
        rmse = testModelRMSE(x,trial,spikecount,posx_bin,model,track_length,dx)
        all_rmse.append(rmse)
        sys.stdout.write('\nfilename-'+str(x)+'-rmse:'+str(rmse))
        sys.stdout.flush()


    # normAllRmse = all_rmse - avgBaselineRMSE


    # Apply a Savitzky-Golay filter to an RMSE line.
    yhat = savgol_filter(all_rmse, 51, 3) # window size 51, polynomial order 3
    plt.clf()
    #plt.plot(plotRange,normAllRmse)
    plt.plot(plotRange,yhat, color='red')
    plt.title('RMSE from Model Trained on Trial 1-50')
    plt.xlabel('Trial')
    plt.ylabel('Normalized RMSE')
    # plt.show()
    plt.savefig('./RMSE/'+filename+'_RMSE.png')
    return {'x':plotRange,'y':yhat, 'model':model}

def circular_rmse(y_test,y_pred,circum):
    error = np.minimum(abs(y_test-y_pred),abs(abs(y_test-y_pred)-circum))
    rmse = math.sqrt(sum(error**2)/len(error))
    return rmse


def testModelRMSE(testTrial, trial,spikecount,posx_bin,model,track_length, dx):
    # test the model 
    test_set = squeeze(trial==testTrial)
    X_test = np.transpose(spikecount[:,test_set])
    y_test = squeeze(posx_bin[test_set])
    y_pred = model.predict(X_test)
    y_pred_proba = model.predict_proba(X_test)

    # compute RMSE
    pred_rmse_bin = circular_rmse(y_test,y_pred,track_length/dx)
    pred_rmse_cm = pred_rmse_bin * dx
    return pred_rmse_cm

##### RUN SCRIPT #####
start_time = time.time()
#seshList = glob.glob('*.mat')
seshList = ['/Users/keimasuda/Desktop/fkm_analysis/combinedSesh/G1_190817_baseline1+controlinjx1+ketamine1.mat']

for sesh in seshList:
    try:
        sessionStr = sesh
        params = {'TrackStart':0,'TrackEnd':400,'TimeBin':0.02}
        logRegDict = trainTestLogisticRegression(sessionStr,params)
        sys.stdout.write('Saving'+sessionStr+'to pickl...\n')
        sys.stdout.flush()
        pickle.dump(logRegDict, open( sessionStr+'_model.pkl', 'wb' ) )
    except Exception as e:
        sys.stdout.write('\nCould not process: '+str(sesh))
        sys.stdout.write(str(e))    
        sys.stdout.flush()

sys.stdout.write("--- %s seconds ---" % (time.time() - start_time))
sys.stdout.flush()

