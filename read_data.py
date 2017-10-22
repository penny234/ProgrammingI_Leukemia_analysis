#this is the script to get the merged_df, training/test sets and cox-score for 21148 genes

import pandas as pd
import numpy as np
import random
from scipy import stats
import matplotlib.pyplot as plt
#1. 
#get the clinical-gene.csv 
cgdata=pd.read_csv("/Users/Penny/Desktop/dataset/ALL_P1/cgdata.csv")
cgdata=cgdata.dropna()  #drop NA reading get 207 patients data
cgdata.shape #get the dimension of the cgdata 207*21154

#2. 
#read the clinical data to see some features of the patients
df=pd.read_excel("/Users/Penny/Desktop/dataset/ALL_P1/clinical.xlsx",sheetname="Clinical Data")
#df contains 255 patients information
some_values=cgdata['patient_id']
sub_df=df.loc[df['TARGET USI'].isin(some_values)] #get subdataframe of df based on 'TARGET USI' ; this mainly contains clinical data
#get the features of clinical data sub_df 
print(sub_df.describe(include='all').transpose()) 
print(sub_df['Race'].value_counts())
print(sub_df['Ethnicity'].value_counts())
print(sub_df['Vital Status'].value_counts())

#3. 
#merged_df, stands for new dataframe which isthe merged dataset we need to use, it contains vital staus 
#
merged_df=pd.merge(cgdata, sub_df[['TARGET USI','Vital Status','Age at Diagnosis in Days','MRD Day 29','WBC at Diagnosis']], left_on='patient_id',right_on='TARGET USI').drop('TARGET USI',1)

#change the 'Age at diagnosis in days' to 'Age at diagnosis in years' 
merged_df['Age at Diagnosis in Days']=merged_df['Age at Diagnosis in Days'].apply(lambda x: x/365)
merged_df=merged_df.rename(columns = {'Age at Diagnosis in Days':'Age at Diagnosis in years'})
merged_df['Age at Diagnosis in years']=[0 if ((x<=10) &(x>=1)) else 1 for x in merged_df['Age at Diagnosis in years']]

#change 'MRD Day 29' to 1 if it is positive; if it is not positive set the value equal 0
merged_df['MRD Day 29']=[1 if x>0 else 0 for x in merged_df['MRD Day 29']] 
#change the 'WBC at Diagnosis', if >5 assign 1 else 0
merged_df['WBC at Diagnosis']=[1 if x>5 else 0 for x in merged_df['WBC at Diagnosis']]


#assign new value 'death' to merged_df
merged_df=merged_df.assign(death=pd.Series(np.ones(len(merged_df['patient_id']))))
merged_df.loc[(merged_df['Vital Status']!='Dead'),'death']=0
merged_df['death']=merged_df['death'].astype(int)


#4. 
#change the gene level V1-V21148 to numeric values 
header=list(merged_df.columns.values) #to see the header of merged_df

#transfer to numeric value 
cols = ['V{0}'.format(element) for element in range(1,21149)]
merged_df.loc[:,cols] = merged_df.loc[:,cols].apply(pd.to_numeric, errors='coerce', axis=1) 




