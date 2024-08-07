#%% Sgro - noise 

import os

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

with open('pop_firingRate_Sgro_noise.txt', 'a') as file:
    for i in range(26):
        for j in range(26):
            line = 'bsub -J "popSgro_n_%.2i_%.2i" -n 2 -W 100 -o /dev/null \'python3 Population_DensityFlowSweep_Sgro_Noise.py %i %i > pop.qlog\'' % (i,j,i,j)
            file.write(line + '\n')
            
#%% Sgro - no noise 

import os

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

with open('pop_firingRate_Sgro_noNoise.txt', 'a') as file:
    for i in range(26):
        for j in range(26):
            line = 'bsub -J "popSgro_n_%.2i_%.2i" -n 2 -W 100 -o /dev/null \'python3 Population_DensityFlowSweep_Sgro_noNoise.py %i %i > pop.qlog\'' % (i,j,i,j)
            file.write(line + '\n')

#%% Goldbeter: noise 

import os

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

with open('pop_firingRate_Goldbeter_noise.txt', 'a') as file:
    for i in range(25):
        for j in range(25):
            line = 'bsub -J "popGB_n_%.2i_%.2i" -n 2 -W 100 -o /dev/null \'python3 Population_DensityFlowSweep_Goldbeter_Noise.py %i %i > pop.qlog\'' % (i,j,i,j)
            file.write(line + '\n')

#%% Goldbeter: no noise 

import os

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

with open('pop_firingRate_Goldbeter_noNoise.txt', 'a') as file:
    for i in range(25):
        for j in range(25):
            line = 'bsub -J "popGB_n_%.2i_%.2i" -n 2 -W 100 -o /dev/null \'python3 Population_DensityFlowSweep_Goldbeter_noNoise.py %i %i > pop.qlog\'' % (i,j,i,j)
            file.write(line + '\n')

#%% Maeda: no noise 

import os

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

with open('pop_firingRate_Maeda_noNoise.txt', 'a') as file:
    for i in range(26):
        for j in range(26):
            line = 'bsub -J "popMaeda_n_%.2i_%.2i" -n 2 -W 100 -o /dev/null \'python3 Population_DensityFlowSweep_Maeda_noNoise.py %i %i > pop.qlog\'' % (i,j,i,j)
            file.write(line + '\n')

#%% Maeda: noise 

import os

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

with open('pop_firingRate_Maeda_Noise.txt', 'a') as file:
    for i in range(26):
        for j in range(26):
            line = 'bsub -J "popMaeda_n_%.2i_%.2i" -n 2 -W 100 -o /dev/null \'python3 Population_DensityFlowSweep_Maeda_Noise.py %i %i > pop.qlog\'' % (i,j,i,j)
            file.write(line + '\n')

#%% Gregor: no noise 

import os

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

with open('pop_firingRate_Gregor_noNoise.txt', 'a') as file:
    for i in range(26):
        for j in range(26):
            line = 'bsub -J "popGregor_n_%.2i_%.2i" -n 2 -W 100 -o /dev/null \'python3 Population_DensityFlowSweep_Gregor_noNoise.py %i %i > pop.qlog\'' % (i,j,i,j)
            file.write(line + '\n')

#%% Gregor: noise 

import os

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

with open('pop_firingRate_Gregor_Noise.txt', 'a') as file:
    for i in range(26):
        for j in range(26):
            line = 'bsub -J "popGregor_n_%.2i_%.2i" -n 2 -W 100 -o /dev/null \'python3 Population_DensityFlowSweep_Gregor_Noise.py %i %i > pop.qlog\'' % (i,j,i,j)
            file.write(line + '\n')

#%% Kamino: no noise 

import os

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

with open('pop_firingRate_Kamino_noNoise.txt', 'a') as file:
    for i in range(31):
        for j in range(31):
            line = 'bsub -J "popKamino_n_%.2i_%.2i" -n 2 -W 100 -o /dev/null \'python3 Population_DensityFlowSweep_Kamino_noNoise.py %i %i > pop.qlog\'' % (i,j,i,j)
            file.write(line + '\n')

#%% Kamino: no noise - log

import os

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

with open('pop_firingRate_Kamino_noNoise_log.txt', 'a') as file:
    for i in range(31):
        for j in range(31):
            line = 'bsub -J "popKamino_n_%.2i_%.2i" -n 2 -W 100 -o /dev/null \'python3 Population_DensityFlowSweep_Kamino_noNoise_log.py %i %i > pop.qlog\'' % (i,j,i,j)
            file.write(line + '\n')
    
#%% Kamino: noise 

import os

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

with open('pop_firingRate_Kamino_Noise.txt', 'a') as file:
    for i in range(31):
        for j in range(31):
            line = 'bsub -J "popKamino_n_%.2i_%.2i" -n 2 -W 100 -o /dev/null \'python3 Population_DensityFlowSweep_Kamino_Noise.py %i %i > pop.qlog\'' % (i,j,i,j)
            file.write(line + '\n')
            
#%% Kamino: noise - log

import os

os.chdir('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2')

with open('pop_firingRate_Kamino_Noise_log.txt', 'a') as file:
    for i in range(31):
        for j in range(31):
            line = 'bsub -J "popKamino_n_%.2i_%.2i" -n 2 -W 100 -o /dev/null \'python3 Population_DensityFlowSweep_Kamino_Noise_log.py %i %i > pop.qlog\'' % (i,j,i,j)
            file.write(line + '\n')