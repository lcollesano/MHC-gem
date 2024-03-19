import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import operator
import random
from scipy.optimize import minimize
import seaborn as sns


aminos=['R','H','K','D','E','S','T', 'N', 'Q', 'C','G', 'P','A','V','I','L', 'M', 'F','Y','W']


#AMINO VECTORS OF 0 AND 1 FOR MATRIX MULTIPLICATION AMINOVECTOR - ENERGY MATRIX

for i in range(len(aminos)):
    k = aminos[i]  
    locals()[str(k)+'_vec']= np.empty(20)
    for j in range(20):
        if j == i:
            locals()[str(k)+'_vec'][j] = 1
        else:
            locals()[str(k)+'_vec'][j] = 0

for i in range(len(aminos)):
    k = aminos[i]  
    locals()[str(k)+'_vec']= np.empty(20)
    for j in range(20):
        if j == i:
            locals()[str(k)+'_vec'][j] = 1
        else:
            locals()[str(k)+'_vec'][j] = 0

def calc_energy(pep):
    all_pep =[]
    for k in range(len(pep)):
        vec =[]
        for i in range(len(pep[k])):
            for j in aminos:
                if pep[k][i] == j:
                    vec.append(globals()[str(j)+'_vec'])
        newpep = (np.array(vec).T)
        #energy.append(sum(en))
        all_pep.append(newpep)
        
    return np.array(all_pep)

def quadratic(pep, Table, Coeff):
    lambd = Coeff
    aa=calc_energy(pep)
    F_lin = []
    for i in range(len(pep)):
        kd = np.sum(aa[i]*Table)
        res= kd + lambd*kd**2
        F_lin.append(float(res))
        
    return F_lin


def linear(pep, Table):
    aa=calc_energy(pep)
    F_lin = []
    for i in range(len(pep)):
        kd = np.sum(aa[i]*Table)
        res= kd 
        F_lin.append(float(res))
        
    return F_lin



########## HERE CHOOSE ALLELES ##########
## FOR MHCI => for HLA-A*1101-> '1101', for HLA-A*0201-> '0201', for HLA-B*0702-> '0702'
## FOR MHCII => for HLA-DRB1*0101-> DRB1010


allele = 'DRB10101'

if allele == '1101':
    lambd= -0.041887053950741604
if allele == '0201':
    lambd= -0.04019313815831509
if allele == '0702':
    lambd= -0.04328069223031524
if allele == 'DRB10101':
    lambd= -0.11692429759043844
###### HERE IMPORT PEPTIDES #####

### only 9mers accepted, in case of MHCII use just cores. 
### Substitute the default file with the list (or pd.DataFrame, or array) with the desired sequences.  

in_file = pd.read_csv('test_out.csv')
sequences = in_file['seq']

### CHOOSE THE MODEL (quadratic or linear)

model = linear

if model == quadratic:
	EM = pd.read_csv('Tables/Global_em_'+str(allele)+'.csv')
	EM = EM.drop(['AA'],axis=1).to_numpy()
	EM = EM.reshape((20,9))
	en = quadratic(sequences, EM, lambd)
	en_res = pd.DataFrame(en)
	en_res.to_csv('out_quadratic.csv', header=False, index=False)

if model == linear:
	EM = pd.read_csv('Tables/Linear_em_'+str(allele)+'.csv')
	EM = EM.drop(['AA'],axis=1).to_numpy()
	EM = EM.reshape((20,9))
	en_linear = linear(sequences, EM)
	en_res = pd.DataFrame(en_linear)
	en_res.to_csv('out_linear.csv', header=False, index=False)







