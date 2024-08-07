# Population Firing Rate Conslidate 

import os
import pickle 
import numpy as np 

#%% Sgro - noise 

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/noise/Sgro')
firingRate_Sgro_noise = np.zeros([26,26])
for file in os.listdir():
    jIndex = int(file.split('_')[1])
    rhoIndex = int(file.split('_')[3].split('.')[0])
    with open(file,'rb') as f:
        ts,As,firingRate,peakHeight = pickle.load(f)
    firingRate_Sgro_noise[rhoIndex,jIndex] = firingRate

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Sgro_noise.pickle','wb') as f:
    pickle.dump(firingRate_Sgro_noise,f,pickle.HIGHEST_PROTOCOL)

#%% Sgro - no noise  

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/no_noise/Sgro')
firingRate_Sgro_noNoise = np.zeros([26,26])
for file in os.listdir():
    jIndex = int(file.split('_')[1])
    rhoIndex = int(file.split('_')[3].split('.')[0])
    with open(file,'rb') as f:
        ts,As,firingRate,peakHeight = pickle.load(f)
    firingRate_Sgro_noNoise[rhoIndex,jIndex] = firingRate

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Sgro_noNoise.pickle','wb') as f:
    pickle.dump(firingRate_Sgro_noNoise,f,pickle.HIGHEST_PROTOCOL)

#%% Goldbeter - noise

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/noise/Goldbeter')
firingRate_Goldbeter_Noise = np.zeros([25,25])
for file in os.listdir():
    kcIndex = int(file.split('_')[1])
    hIndex = int(file.split('_')[3].split('.')[0])
    with open(file,'rb') as f:
        ts,bs,firingRate,peakHeight = pickle.load(f)
    firingRate_Goldbeter_Noise[hIndex,kcIndex] = firingRate

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Goldbeter_Noise.pickle','wb') as f:
    pickle.dump(firingRate_Goldbeter_Noise,f,pickle.HIGHEST_PROTOCOL)

#%% Goldbeter - no noise

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/no_noise/Goldbeter')
firingRate_Goldbeter_noNoise = np.zeros([25,25])
for file in os.listdir():
    kcIndex = int(file.split('_')[1])
    hIndex = int(file.split('_')[3].split('.')[0])
    with open(file,'rb') as f:
        ts,bs,firingRate,peakHeight = pickle.load(f)
    firingRate_Goldbeter_noNoise[hIndex,kcIndex] = firingRate

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Goldbeter_noNoise.pickle','wb') as f:
    pickle.dump(firingRate_Goldbeter_noNoise,f,pickle.HIGHEST_PROTOCOL)
    
#%% Maeda - no noise

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/no_noise/Maeda')
firingRate_Maeda_noNoise = np.zeros([26,26])
for file in os.listdir():
    gammaIndex = int(file.split('_')[1])
    rhoIndex = int(file.split('_')[3].split('.')[0])
    with open(file,'rb') as f:
        ts,cAMPis,firingRate,peakHeight = pickle.load(f)
    firingRate_Maeda_noNoise[rhoIndex,gammaIndex] = firingRate

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Maeda_noNoise.pickle','wb') as f:
    pickle.dump(firingRate_Maeda_noNoise,f,pickle.HIGHEST_PROTOCOL)
    
#%% Maeda - noise

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/noise/Maeda')
firingRate_Maeda_Noise = np.zeros([26,26])
for file in os.listdir():
    gammaIndex = int(file.split('_')[1])
    rhoIndex = int(file.split('_')[3].split('.')[0])
    with open(file,'rb') as f:
        ts,cAMPis,firingRate,peakHeight = pickle.load(f)
    firingRate_Maeda_Noise[rhoIndex,gammaIndex] = firingRate

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Maeda_Noise.pickle','wb') as f:
    pickle.dump(firingRate_Maeda_Noise,f,pickle.HIGHEST_PROTOCOL)

#%% Gregor - no noise

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/no_noise/Gregor')
firingRate_Gregor_noNoise = np.zeros([26,26])
for file in os.listdir():
    rhoIndex = int(file.split('_')[1])
    kIndex = int(file.split('_')[3].split('.')[0])
    with open(file,'rb') as f:
        ts,cAMPis,firingRate,peakHeight = pickle.load(f)
    firingRate_Gregor_noNoise[rhoIndex,kIndex] = firingRate

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Gregor_noNoise.pickle','wb') as f:
    pickle.dump(firingRate_Gregor_noNoise,f,pickle.HIGHEST_PROTOCOL)
    
#%% Gregor - noise

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/noise/Gregor')
firingRate_Gregor_Noise = np.zeros([26,26])
for file in os.listdir():
    rhoIndex = int(file.split('_')[1])
    kIndex = int(file.split('_')[3].split('.')[0])
    with open(file,'rb') as f:
        ts,cAMPis,firingRate,peakHeight = pickle.load(f)
    firingRate_Gregor_Noise[rhoIndex,kIndex] = firingRate

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Gregor_Noise.pickle','wb') as f:
    pickle.dump(firingRate_Gregor_Noise,f,pickle.HIGHEST_PROTOCOL)
    
#%% Kamino - no noise

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/no_noise/Kamino')
firingRate_Kamino_noNoise = np.zeros([31,31])
for file in os.listdir():
    gammaIndex = int(file.split('_')[1])
    rhoIndex = int(file.split('_')[3].split('.')[0])
    with open(file,'rb') as f:
        ts,cAMPis,firingRate,peakHeight = pickle.load(f)
    firingRate_Kamino_noNoise[rhoIndex,gammaIndex] = firingRate

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Kamino_noNoise.pickle','wb') as f:
    pickle.dump(firingRate_Kamino_noNoise,f,pickle.HIGHEST_PROTOCOL)
    
#%% Kamino - no noise, log

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/no_noise/Kamino_log')
rates = np.zeros([31,31])
for file in os.listdir():
    gammaIndex = int(file.split('_')[1])
    rhoIndex = int(file.split('_')[3].split('.')[0])
    with open(file,'rb') as f:
        ts,cAMPis,firingRate,peakHeight = pickle.load(f)
    rates[rhoIndex,gammaIndex] = firingRate

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Kamino_noNoise_log.pickle','wb') as f:
    pickle.dump(rates,f,pickle.HIGHEST_PROTOCOL)
    
#%% Kamino - noise

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/noise/Kamino')
rates = np.zeros([31,31])
for file in os.listdir():
    gammaIndex = int(file.split('_')[1])
    rhoIndex = int(file.split('_')[3].split('.')[0])
    with open(file,'rb') as f:
        ts,cAMPis,firingRate,peakHeight = pickle.load(f)
    rates[rhoIndex,gammaIndex] = firingRate

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Kamino_Noise.pickle','wb') as f:
    pickle.dump(rates,f,pickle.HIGHEST_PROTOCOL)
    
#%% Kamino - noise, log

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/noise/Kamino_log')
rates = np.zeros([31,31])
for file in os.listdir():
    gammaIndex = int(file.split('_')[1])
    rhoIndex = int(file.split('_')[3].split('.')[0])
    with open(file,'rb') as f:
        ts,cAMPis,firingRate,peakHeight = pickle.load(f)
    rates[rhoIndex,gammaIndex] = firingRate

with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/pop_firingRate_data/firingRates_Kamino_Noise_log.pickle','wb') as f:
    pickle.dump(rates,f,pickle.HIGHEST_PROTOCOL)