import os,re,math
import pandas as pd
import numpy as np
from scipy.spatial import KDTree
from math import factorial as ft
from math import log

class dump:
    def __init__(self, time_step, xyz_df):
        self.time_step = time_step
        self.len = len(xyz_df)
        self.xyz_df = xyz_df


    def get_xyz_type(self, type_num):
        df_selected = self.xyz_df[self.xyz_df[1] == type_num]
        return df_selected

def read_dumps(file_path):
    line_pointer = 0
    num_of_atom_df = pd.read_csv(file_path, skiprows=line_pointer + 3, \
                            nrows=1, header=None)
    num_of_atom = int(num_of_atom_df[0])

    dump_list = []
    while True:
        try:
            time_step_df = pd.read_csv(file_path, skiprows=line_pointer + 1, \
                                    nrows=1, header=None)
        except:
            break
        time_step = int(time_step_df[0])
        #print(time_step)
        xyz_df = pd.read_csv(file_path, skiprows=line_pointer + 9, \
                                   nrows=num_of_atom, sep=' ', header=None)
        #print(xyz_df)
#------------------------- total frames in trj / skip = how many frames read in -----------------------------------#
        skip = 1
        print('Time Step: ', time_step, ' completed.')
        line_pointer += (9 + num_of_atom) * skip


        dump_list.append(dump(time_step, xyz_df))
    return dump_list

file_path='/data/home/djy5/Work/HHe/solar/Jupiter/P2T12/'
name='P2T12'
nNB=200

dump_list= read_dumps(file_path=file_path+name+'.lammpstrj')
system_box_df = pd.read_csv(file_path+name+'.lammpstrj', skiprows= 5,nrows=3, \
                            header=None, delim_whitespace=True).values.astype(np.float)

xl = system_box_df[0][1] - system_box_df[0][0]

#--------------------------------------grow cubic-----------------------------------------#

all_frame = len(dump_list)
Natoms=dump_list[0].xyz_df.shape[0]
frame=0
time_list = []
all_sHe = []
all_sH = []
while frame < all_frame:
    print('Now the calculating frame is:' + str(frame))
    current_frame=dump_list[frame]
    TYPE = current_frame.xyz_df[[1]]
    XYZ = current_frame.xyz_df[[2,3,4]]
    xyz=np.array(np.array(XYZ).tolist()) #free XYZ from dataframe to an array of sublists
    kdtree=KDTree(xyz)
    dist,points=kdtree.query(xyz,nNB)
    df = pd.DataFrame(columns=['nbH','nbHe','xiH','xiHe'])
    Type = np.array(np.array(TYPE).tolist())
    np.reshape(Type[points], (nNB, Natoms))
    xi = pd.DataFrame(pd.DataFrame(np.reshape(Type[points], (Natoms, nNB)))[0])
    df_nb = pd.DataFrame(np.reshape(Type[points], (Natoms, nNB))).drop([0], axis=1)
    xi['nbH'] = (df_nb == 1).sum(axis=1)
    xi['nbHe'] = (df_nb == 2).sum(axis=1)
    xi_1 = list(xi.groupby(0))[0][1].groupby(['nbH', 'nbHe']).size()
    xi_2 = list(xi.groupby(0))[1][1].groupby(['nbH', 'nbHe']).size()
    xi_all = xi.groupby(['nbH', 'nbHe']).size()
    all_sHe.append(np.array(np.repeat((xi_1 / xi_all).dropna(), xi_1)))
    all_sH.append(np.array(np.repeat((xi_2 / xi_all).dropna(), xi_2)))
    time_list.append(frame * 2)  # unit:ps, postion dump every 10000 steps,dt=0.2fs*10000=2ps then * skip=2ps
    frame = frame + 1
sHe = pd.DataFrame(all_sHe).T
sHe.columns = time_list
sH = pd.DataFrame(all_sH).T
sH.columns = time_list
sH.to_csv('pH_' + name +  '.csv')
sHe.to_csv('pHe_' + name +  '.csv')
