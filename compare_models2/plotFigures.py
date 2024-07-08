# plot figures for Chuqiao's paper 

import pickle 
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

#%% Figure 3 (population add cAMP): pull data

# experimental data
my_dir = '//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/exp_data/'
Sgro2015Figure6excel = pd.read_excel(my_dir+r'Sgro2015DataFormattedforPython.xlsx',sheet_name='Figure6')

# Sgro 
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/popAddcAMP_data/popAddcAMP_Sgro.pickle','rb') as f:
    ts_Sgro,As_Sgro,As_sc_Sgro = pickle.load(f)

# Gregor
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/popAddcAMP_data/popAddcAMP_Gregor.pickle','rb') as f:
    ts,thetais,cAMPes,cAMPis_sc,cAMPis = pickle.load(f)

# Goldbeter
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/popAddcAMP_data/popAddcAMP_Goldbeter.pickle','rb') as f:
    ts,ps,bs,gs = pickle.load(f)
    
# Maeda
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/popAddcAMP_data/popAddcAMP_Maeda.pickle','rb') as f:
    t,cAMPi_traces = pickle.load(f)

# Kamino
with open('//prfs.hhmi.org/sgrolab/mark/dicty_proj/dictymodels/compare_models2/popAddcAMP_data/popAddcAMP_Kamino.pickle','rb') as f:
    t,x_traces,y_traces,z_traces = pickle.load(f)


#%% plot 

mycolors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']     
          
# label_font_size = 18; trace_width = 3; tick_font_size=15

abcd_font_size = 28
label_font_size=24
title_font_size = 26
sublabel_font_size = 22
trace_width=3
tick_font_size=20

fig5 = plt.figure(figsize=(11, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.3)

#void384 data = MATLAB structure, need to figure out what is what

#ax01= fig5.add_subplot(grid[0, 0])
#ax01.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["No External cAMP Mean Trace"],
#                               color = 'k', linewidth=trace_width)
#ax01.axvspan(60, 120, alpha=0.2, color='b')
#ax01.set_ylim([-0.1,0.6]);ax01.set_xlim([0,120])
#ax01.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
#ax01.text(0.7,0.85,' No External cAMP', horizontalalignment='center',verticalalignment='center',
#         transform = ax01.transAxes, color = 'k', fontsize=tick_font_size)

ax02= fig5.add_subplot(grid[0, 0])
ax02.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["Low External cAMP Mean Trace"],
                               color = 'k', linewidth=trace_width)
ax02.axvspan(60, 120, alpha=0.2, color='b')
ax02.set_ylim([-0.1,0.6]);ax02.set_xlim([0,120])
ax02.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax02.text(0.7,0.75,' Low External cAMP, \n 5-10nM', horizontalalignment='center',verticalalignment='center',
         transform = ax02.transAxes, color = 'k', fontsize=tick_font_size)
ax02.set_title('Experiment',size=title_font_size)
ax03= fig5.add_subplot(grid[1, 0])
ax03.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["Intermediate External cAMP Mean Trace"],
                               color = 'k', linewidth=trace_width)
ax03.axvspan(60, 120, alpha=0.2, color='b')
ax03.set_ylim([-0.1,0.6]);ax03.set_xlim([0,120])
ax03.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax03.text(0.7,0.75,' Intermediate External \n cAMP, 10-20nM', horizontalalignment='center',verticalalignment='center',
         transform = ax03.transAxes, color = 'k', fontsize=tick_font_size)

ax04= fig5.add_subplot(grid[2, 0])
ax04.plot(Sgro2015Figure6excel["Times (min)"],Sgro2015Figure6excel["High External cAMP Mean Trace"],
                               color = 'k', linewidth=trace_width)
ax04.axvspan(60, 120, alpha=0.2, color='b')
ax04.set_ylim([-0.1,0.6]);ax04.set_xlim([0,120])
ax04.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax04.text(0.7,0.75,' High External cAMP, \n 100nM', horizontalalignment='center',verticalalignment='center',
         transform = ax04.transAxes, color = 'k', fontsize=tick_font_size)

fig5.text(0.02, 0.9, 'A', color='b', fontsize=abcd_font_size, ha='center')

fig5.text(0.02, 0.5, r'FRET Signal, A.U.', fontsize=label_font_size,va='center', rotation='vertical')
fig5.text(0.5, 0.05, 'Time (min)', fontsize=label_font_size, ha='center')


#% Plot  3 traces: low, medium and high [cAMP]ext 

fig = plt.figure(figsize=(11, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.3)

ax1= fig.add_subplot(grid[0, 0])
for count in range(5):
    ax1.plot(t_plot_Sgro,A_traces_single_cell[0,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax1.plot(t_plot_Sgro,A_traces[0,:], color=mycolors[5],linewidth=trace_width)
ax1.text(0.7,0.9,r'Low $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax1.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax1.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax1.axvspan(15, 30, alpha=0.2, color='g')
ax1.set_xlim([0,30]); ax1.set_ylim([-0.25,1.25])

ax2= fig.add_subplot(grid[1, 0])
for count in range(5):
    ax2.plot(t_plot_Sgro,A_traces_single_cell[1,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax2.plot(t_plot_Sgro,A_traces[1,:], color=mycolors[5],linewidth=trace_width)
    

ax2.text(0.7,0.9,r'Intermediate $cAMP_{e}$'+' input', horizontalalignment='center',verticalalignment='center',
     transform = ax2.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax2.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax2.axvspan(15, 30, alpha=0.2, color='g')
ax2.set_xlim([0,30]); ax2.set_ylim([-0.25,1.25])

ax3= fig.add_subplot(grid[2, 0])
for count in range(5):
    ax3.plot(t_plot_Sgro,A_traces_single_cell[2,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax3.plot(t_plot_Sgro,A_traces[2,:], color=mycolors[5],linewidth=trace_width)
ax3.text(0.7,0.9,r'High $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax3.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax3.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax3.axvspan(15, 30, alpha=0.2, color='g')
ax3.set_xlim([0,30]); ax3.set_ylim([-0.25,1.25])

fig.text(0.02, 0.9, 'E', color='g', fontsize=abcd_font_size, ha='center')
fig.text(0.5, 0.04, 'Time, A.U.',fontsize=label_font_size, ha='center')
fig.text(0.02, 0.5, r'$cAMP_{i}$',fontsize=label_font_size, va='center', rotation='vertical')
fig.text(0.5, 0.9, 'Sgro 2015',color = mycolors[5],fontsize=title_font_size, ha='center')
plt.show()

# Gregor 
#% Plot  3 traces: low, medium and high [cAMP]ext 
fig = plt.figure(figsize=(11, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.3)

ax1= fig.add_subplot(grid[0, 0])
for count in range(10):
    campCyto_traces_single_cell[0,count,:] = campCyto_traces_single_cell[0,count,:]# /np.amax(campCyto_traces_single_cell[0,count,:])
    ax1.plot(t_plot_Gregor,campCyto_traces_single_cell[0,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax1.plot(t_plot_Gregor,campCyto_traces[0,:], alpha=0.8, color=mycolors[2],linewidth=trace_width+1)  

ax1.text(0.7,0.9,r'Low $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax1.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax1.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax1.axvspan(15, 30, alpha=0.2, color='g')
ax1.set_xlim([t_plot_Gregor[1],30]); ax1.set_ylim([-0.1,1.3])

ax2= fig.add_subplot(grid[1, 0])
for count in range(10):
    campCyto_traces_single_cell[1,count,:] = campCyto_traces_single_cell[1,count,:]# /np.amax(campCyto_traces_single_cell[1,count,:])
    ax2.plot(t_plot_Gregor,campCyto_traces_single_cell[1,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax2.plot(t_plot_Gregor,campCyto_traces[1,:], alpha=0.8, color=mycolors[2],linewidth=trace_width+1)

ax2.text(0.7,0.9,r'Intermediate $cAMP_{e}$'+' input', horizontalalignment='center',verticalalignment='center',
     transform = ax2.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax2.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax2.axvspan(15, 30, alpha=0.2, color='g')
ax2.set_xlim([t_plot_Gregor[1],30]); ax2.set_ylim([-0.1,1.3])

ax3= fig.add_subplot(grid[2, 0])
for count in range(10):
    campCyto_traces_single_cell[2,count,:] = campCyto_traces_single_cell[2,count,:]# /np.amax(campCyto_traces_single_cell[2,count,:])
    ax3.plot(t_plot_Gregor,campCyto_traces_single_cell[2,count,:],
             color='darkgrey',alpha=0.5, linewidth=trace_width-1)
ax3.plot(t_plot_Gregor,campCyto_traces[2,:], color=mycolors[2],alpha = 0.8, linewidth=trace_width+1)

ax3.text(0.7,0.9,r'High $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax3.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax3.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax3.axvspan(15, 30, alpha=0.2, color='g')
ax3.set_xlim([t_plot_Gregor[1],30]); ax3.set_ylim([-0.1,1.3])

fig.text(0.02, 0.9, 'D', color='g', fontsize=abcd_font_size, ha='center')
fig.text(0.5, 0.04, 'Time, A.U.',fontsize=label_font_size, ha='center')
fig.text(0.02, 0.5, r'$cAMP_{i}$',fontsize=label_font_size, va='center', rotation='vertical')
fig.text(0.5, 0.9, 'Gregor 2010',color = mycolors[2],fontsize=label_font_size, ha='center')
plt.show()

# GOLDBETER

fig = plt.figure(figsize=(11, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.3)

ax1= fig.add_subplot(grid[0, 0])
ax1.plot(t_plot_Goldbeter,b_traces[0,:], color=mycolors[0],linewidth=trace_width)
ax1.text(0.7,0.9,r'Low $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax1.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax1.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax1.axvspan(15, 30, alpha=0.2, color='g')
ax1.set_xlim([0,30]); ax1.set_ylim([-0.25,1.5])

ax2= fig.add_subplot(grid[1, 0])
ax2.plot(t_plot_Goldbeter,b_traces[1,:], color=mycolors[0],linewidth=trace_width)
ax2.text(0.7,0.8,r'Intermediate $cAMP_{e}$'+'\n input', horizontalalignment='center',verticalalignment='center',
     transform = ax2.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax2.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax2.axvspan(15, 30, alpha=0.2, color='g')
ax2.set_xlim([0,30]); ax2.set_ylim([-0.25,1.5])

ax3= fig.add_subplot(grid[2, 0])
ax3.plot(t_plot_Goldbeter,b_traces[2,:], color=mycolors[0],linewidth=trace_width)
ax3.text(0.7,0.8,r'High $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax3.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax3.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax3.axvspan(15, 30, alpha=0.2, color='g')
ax3.set_xlim([0,30]); ax3.set_ylim([-0.25,1.5])

fig.text(0.02, 0.9, 'B', color='g', fontsize=abcd_font_size, ha='center')
fig.text(0.5, 0.04, 'Time, A.U.',fontsize=label_font_size, ha='center')
fig.text(0.02, 0.5, r'$cAMP_{i}$',fontsize=label_font_size, va='center', rotation='vertical')
fig.text(0.5, 0.9, 'Martiel 1987',color = mycolors[0],fontsize=label_font_size, ha='center')
plt.show()


# MAEDA

fig = plt.figure(figsize=(11, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.3)

ax1= fig.add_subplot(grid[0, 0])
ax1.plot(t_plot_Maeda,cAMPi_traces[0,:], color=mycolors[1],linewidth=trace_width)
ax1.text(0.7,0.9,r'Low $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax1.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax1.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax1.axvspan(15, 30, alpha=0.2, color='g')
ax1.set_xlim([0,30]); ax1.set_ylim([0.1,0.9])

ax2= fig.add_subplot(grid[1, 0])
ax2.plot(t_plot_Maeda, cAMPi_traces[1,:], color=mycolors[1],linewidth=trace_width)
ax2.text(0.7,0.8,r'Intermediate $cAMP_{e}$'+'\n input', horizontalalignment='center',verticalalignment='center',
     transform = ax2.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax2.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax2.axvspan(15, 30, alpha=0.2, color='g')
ax2.set_xlim([0,30]); ax2.set_ylim([0.1,0.9])

ax3= fig.add_subplot(grid[2, 0])
ax3.plot(t_plot_Maeda, cAMPi_traces[2,:], color=mycolors[1],linewidth=trace_width)
ax3.text(0.7,0.8,r'High $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax3.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax3.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax3.axvspan(15, 30, alpha=0.2, color='g')
ax3.set_xlim([0,30]); ax3.set_ylim([0.1,0.9])

fig.text(0.02, 0.9, 'C', color='g', fontsize=abcd_font_size, ha='center')
fig.text(0.5, 0.04, 'Time, A.U.',fontsize=label_font_size, ha='center')
fig.text(0.02, 0.5, r'$cAMP_{i}$',fontsize=label_font_size, va='center', rotation='vertical')
fig.text(0.5, 0.9, 'Maeda 2004',color = mycolors[1],fontsize=label_font_size, ha='center')
plt.show()

# KAMINO

fig = plt.figure(figsize=(11, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.3)

ax1= fig.add_subplot(grid[0, 0])
ax1.plot(t_plot_Kamino, y_traces[0,:], color=mycolors[7],linewidth=trace_width)
ax1.text(0.7,0.9,r'Low $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax1.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax1.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax1.axvspan(15, 30, alpha=0.2, color='g')
ax1.set_xlim([0,30]); ax1.set_ylim([-0.2,1.2])

ax2= fig.add_subplot(grid[1, 0])
ax2.plot(t_plot_Kamino, y_traces[1,:], color=mycolors[7],linewidth=trace_width)
ax2.text(0.7,0.8,r'Intermediate $cAMP_{e}$'+'\n input', horizontalalignment='center',verticalalignment='center',
     transform = ax2.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax2.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax2.axvspan(15, 30, alpha=0.2, color='g')
ax2.set_xlim([0,30]); ax2.set_ylim([-0.2,1.2])

ax3= fig.add_subplot(grid[2, 0])
ax3.plot(t_plot_Kamino, y_traces[2,:], color=mycolors[7],linewidth=trace_width)
ax3.text(0.7,0.8,r'High $cAMP_{e}$ input', horizontalalignment='center',verticalalignment='center',
     transform = ax3.transAxes, color = 'k', fontsize=tick_font_size)
#    ax.set_xlabel(r'$cAMP_{ext}$ input='+str(alphafval_arr[i])+ 'nM', fontsize=label_font_size)
ax3.tick_params(grid_linewidth = 15, labelsize = tick_font_size)
ax3.axvspan(15, 30, alpha=0.2, color='g')
ax3.set_xlim([0,30]); ax3.set_ylim([-0.2,1.2])

fig.text(0.02, 0.9, 'F', color='g', fontsize=abcd_font_size, ha='center')
fig.text(0.5, 0.04, 'Time, A.U.',fontsize=label_font_size, ha='center')
fig.text(0.02, 0.5, r'$cAMP_{i}$',fontsize=label_font_size, va='center', rotation='vertical')
fig.text(0.5, 0.9, 'Kamino 2017',color = mycolors[7],fontsize=label_font_size, ha='center')
plt.show()