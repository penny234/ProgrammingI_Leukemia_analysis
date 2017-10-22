#this script consider the pca with standardized gene value

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from sklearn import decomposition
from sklearn.decomposition import PCA as sklearnPCA
from lifelines import CoxPHFitter
from getcox import h
from sklearn import preprocessing



def f(train,threshold,test):
    hi=h(train)
    h_score=pd.DataFrame(hi, index=np.array(range(1,21149)))
    gene_ls=h_score.index[h_score.iloc[:,0]>1].tolist()
    candidate_genes=['V{0}'.format(element) for element in gene_ls]

    # qualified genes were selected 

    stdsc = preprocessing.StandardScaler()
    np_scaled_train = stdsc.fit_transform(train.loc[:,candidate_genes])
    np_scaled_test  = stdsc.transform(test.loc[:,candidate_genes])
    pca = sklearnPCA(n_components=1)   
    X_train_pca = pca.fit_transform(np_scaled_train) # This is the result 
    X_test_pca  = pca.transform(np_scaled_test)
    eigen_val=pca.explained_variance_  #eigen value is the explained variance 

    
    #assign pca score to the test dataset 
    test=test.assign(w=pd.Series(np.ones(len(test.patient_id))))
    test['w']=X_test_pca
    testset_surv=test[['event_free_survival_time_days','death','w']]
    
    #do cox-regression

    # Using Cox Proportional Hazards model
    cph = CoxPHFitter()
    cph.fit(testset_surv,'event_free_survival_time_days',event_col='death')
    
    return cph.print_summary()