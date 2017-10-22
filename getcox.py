#latest working version 
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def h(inputdata):
    
#1. get inputdata
#2. get event time
    event_time=inputdata[inputdata['death']==1]['event_free_survival_time_days'].sort_values().unique()
    
    cols = ['V{0}'.format(element) for element in range(1,21149)] #has 21148 genes generate the genenames

    ls=[]
    for i in event_time:
        m=(inputdata['event_free_survival_time_days']>i).sum()
        d=(inputdata['event_free_survival_time_days']==i).sum()
        who=inputdata.index[inputdata['event_free_survival_time_days']==i].tolist()
        ls.append([m,d,who])
    
#the index of at risk set 
    indx=[[] for _ in range(len(event_time))]
    for i in range(len(event_time)):
        indx[i]=inputdata[inputdata['event_free_survival_time_days']>event_time[i]].index.values.tolist()

    
#get Xik
    sum_gene_val=pd.DataFrame(0, index=cols, columns=np.array(range(len(event_time))))
    for i in range(len(event_time)): 
        sum_gene_val.loc[:,i]=inputdata.loc[indx[i],cols].sum()

#get Xik_bar
    ser1=[item[0] for item in ls]
    func1 = lambda x: np.asarray(x) / np.asarray(ser1)
    avg_gene_val=sum_gene_val.apply(func1, axis=1)


#get X*ik
    subj_time=[item[2] for item in ls]
    subj_gene_val=pd.DataFrame(0, index=cols, columns=np.array(range(len(event_time))))
    for i in range(len(event_time)):    
        subj_gene_val.loc[:,i]=inputdata.loc[subj_time[i][0],cols]  #*bug
    
#get dk*Xik_bar
    ser2=[item[1] for item in ls]
    func2 = lambda x: np.asarray(x)*np.asarray(ser2)
    davg_gene_val=avg_gene_val.apply(func2, axis=1)



#get ri
    r=davg_gene_val.subtract(subj_gene_val, axis='columns').sum(axis=1)

    #get Si
    S_gene_val=pd.DataFrame(0, index=cols, columns=np.array(range(len(event_time))))
    temp=inputdata.loc[:,cols].T
    for i in range(len(event_time)):
        S_gene_val.loc[:,i]=(temp.loc[:,indx[i]].sub(avg_gene_val.loc[:,i],axis=0)**2).sum(axis=1)*np.asarray(ser2[i])/np.asarray(ser1[i])

    
    si=S_gene_val.sum(axis=1)**(1/2)
    s0=si.median()    
    hi=np.asarray(r)/(np.asarray(si)+s0)
    hi=np.absolute(hi)  #get absolute value of cox-score

    stats.describe(hi) #mean=1.1501915210128424, variance=0.2270225498339461,  minmax=(2.8752108009228312e-06, 4.0052250126097571)
    np.histogram(hi)   
    plt.hist(hi,bins=20)
    plt.title("Histogram of absolute value of cox-score")
    plt.show()

    return hi







