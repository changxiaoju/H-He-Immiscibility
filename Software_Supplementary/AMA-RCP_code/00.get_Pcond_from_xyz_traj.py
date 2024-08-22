import os,re,math,glob,sys,time
import pandas as pd
import numpy as np
from scipy.spatial import KDTree
from math import factorial as ft
from math import log
sys.path.append('/ADD_PATH_TO/MolecularModTools/OutputInfo/') #on https://github.com/changxiaoju/MolecularModTools/releases/tag/v1.0.0
from ReadBox import read_xyz_frame

file_path='/PATH_TO_XYZ_FILE/'
name='NAME_OF_XYZ_FILE'
nNB=100 #the number of nearest atoms for a central atom

frames = []
with open(file_path+name+'.xyz', 'r')  as xyz_file:
    while True:
        try:
            atom_types,coordinates,cell_length= read_xyz_frame(xyz_file)
            frames.append((coordinates,atom_types))
        except ValueError:
            break

Natoms=len(coordinates)
time_list = []
all_Pcond_A,all_Pcond_B = [],[]
for index,frame in enumerate(frames):
    print('Now the calculating frame is:' + str(index))
    current_frame=frame[0]
    XYZ = current_frame
    xyz=np.array(np.array(XYZ).tolist())
    kdtree=KDTree(xyz)
    dist,points=kdtree.query(xyz,nNB)
    df = pd.DataFrame(columns=['nbA','nbB','xiA','xiB'])
    Type = np.array(frame[1])
    
    #reshape data to Natoms*nNB dataframe
    df_type = pd.DataFrame(np.reshape(Type[points], (Natoms, nNB))).astype('int')
    xi = pd.DataFrame(df_type[0])
    df_nb = df_type.drop([0], axis=1)

    #count the number of A and B atoms in the neigborhood of the ith atom
    xi['nbA'] = (df_nb == 1).sum(axis=1)
    xi['nbB'] = (df_nb == 2).sum(axis=1)

    xi_A = list(xi.groupby(0))[0][1].groupby(['nbA', 'nbB']).size()
    xi_B = list(xi.groupby(0))[1][1].groupby(['nbA', 'nbB']).size()
    xi_all = xi.groupby(['nbA', 'nbB']).size()
    
    #calculate every atom's Pcond
    all_Pcond_A.append(np.array(np.repeat((xi_A / xi_all).dropna(), xi_A)))
    all_Pcond_B.append(np.array(np.repeat((xi_B / xi_all).dropna(), xi_B)))
    
    time_list.append(index * 1) # unit:ps, postion dump every 5000 steps,dt=0.1fs*5000=1ps

Pcond_A = pd.DataFrame(all_Pcond_A).T
Pcond_A.columns = time_list
Pcond_B = pd.DataFrame(all_Pcond_B).T
Pcond_B.columns = time_list
Pcond_A.to_csv('pA_' + name +  '.csv')
Pcond_B.to_csv('pB_' + name +  '.csv')
