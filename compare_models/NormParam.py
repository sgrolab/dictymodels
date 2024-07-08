# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 11:00:54 2020

@author: Chuqiao Huyan

Time and cAMP maginitude normalization parameters. These parameters are kept
constant across the whole paper for each model.

"""
import pickle

Nt_Gregor = 6 
Nt_Sgro = 27
Nt_Goldbeter = 6.94
Nt_Maeda = 3.57
Nt_Kamino = 5.23

Nh_Gregor = 19.6
Nh_Sgro = 3.5
Nh_Goldbeter = 210.53
Nh_Maeda = 3.15 
Nh_Kamino = 0.26

Nh_Kamino_offset = 0.0585
Nh_Sgro_offset = -1.5

with open('normParams.pickle','wb') as f:
    pickle.dump([Nt_Gregor,Nt_Sgro,Nt_Goldbeter,Nt_Maeda,Nt_Kamino,Nh_Gregor,Nh_Sgro,Nh_Goldbeter,Nh_Maeda,Nh_Kamino,Nh_Kamino_offset,Nh_Sgro_offset],f,pickle.HIGHEST_PROTOCOL)