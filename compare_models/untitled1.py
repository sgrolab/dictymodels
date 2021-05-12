# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 12:29:58 2020

@author: ellin
"""



#%%
title_font_size = 20
label_font_size = 20
sublabel_font_size = 18
tick_font_size = 18
trace_width = 2
abcd_font_size = 28

fig = plt.figure(figsize=(cm2inch(48),cm2inch(24)))
grid = plt.GridSpec(3,5, wspace=0.35, hspace=0.6)

axB0= fig.add_subplot(grid[0,0])
axB0.set_xticks([0,2,4,6,8]); 
axB0.set_xticklabels([1,3,5,7,9],fontsize=tick_font_size)
# axB0.set_yticks([0,1,2,3,4,5,6,7]); 
# axB0.set_yticklabels(['1/2','1/4','1/8','1/16','1/32','1/64','1/128'],fontsize=tick_font_size-3)
axB0.set_yticks([0,2,4,6]); 
axB0.set_yticklabels(['1/2','1/8','1/32','1/128'],fontsize=tick_font_size)

axB0.set_title('Experiment', fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
axB0.set_xlabel('Flow Rate (mL/min)', size=sublabel_font_size)
axB0.set_ylabel('Cell Density(mML)', size=sublabel_font_size)
heatmap = axB0.imshow(PopRateExp, cmap='jet') # cmap='jet'
x=[3.5,4.5,5.5,7.5,9.5]
[axB0.axvline(_x, color='white',linewidth=trace_width) for _x in x]
# heatmap.set_clim(0,0.16)
cbar=fig.colorbar(heatmap, ax=axB0,ticks=[0,0.05,0.1,0.15]);
cbar.ax.tick_params(labelsize = tick_font_size) 
cbar.set_label( 'cAMP pulses/min',size=tick_font_size)
#axB0.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB0.text(-0.6 , 1.6, 'a',
         ha='center',va='center',
         transform = axB0.transAxes, color = expcolor, fontsize=abcd_font_size)

# axB1= fig.add_subplot(grid[1,0],xticks=[0,25,50,75,100])
# heatmap = axB1.pcolor(kc_arr_Gold_noise10, oneoverh_arr_Gold_noise10, pop_rate_Gold_noise10.transpose(), cmap='jet') # cmap='jet'
# # axB1.set_xscale('log');
# axB1.set_yscale('log')
# heatmap.set_clim(0,1.5)
# cbar=fig.colorbar(heatmap, ax=axB1, ticks=[0,0.5,1,1.5])
# cbar.ax.tick_params(labelsize = tick_font_size) 
# axB1.tick_params(axis='both', which='major', labelsize=tick_font_size)
# axB1.set_title('Martiel 1987\nwith noise', color = goldcolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# axB1.text(-0.6 , 1.2, 'B',
#          ha='center',va='center',
#          transform = axB1.transAxes, color = simcolor, fontsize=abcd_font_size)

axB1lower= fig.add_subplot(grid[2,0],xticks=[0,25,50,75,100])
heatmap = axB1lower.pcolor(kc_arr_Gold, oneoverh_arr_Gold, pop_rate_Gold.transpose(), cmap='jet') # cmap='jet'
# axB1.set_xscale('log');
axB1lower.set_yscale('log')
heatmap.set_clim(0,1.5)
cbar=fig.colorbar(heatmap, ax=axB1lower, ticks=[0,0.5,1,1.5])
cbar.ax.tick_params(labelsize = tick_font_size) 
axB1lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB1lower.set_title('Martiel 1987\nw/o noise', color = goldcolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

axB1lower.text(-0.6 , 1.2, '(b)',
          ha='center',va='center',
          transform = axB1lower.transAxes, color = simcolor, fontsize=abcd_font_size)

# axB2= fig.add_subplot(grid[1,1],xticks=[0,25,50,75,100])
# heatmap = axB2.pcolor(gamma_arr_Maeda_noise1, rho_arr_Maeda_noise1, pop_rate_Maeda_noise1.transpose(), cmap='jet') # cmap='jet'
# heatmap.set_clim(0,0.55)
# axB2.set_yscale('log')
# cbar=fig.colorbar(heatmap, ax=axB2);cbar.ax.tick_params(labelsize = tick_font_size) 
# axB2.tick_params(axis='both', which='major', labelsize=tick_font_size)
# axB2.set_title('Maeda 2004\nwith noise', color = maedacolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# axB2.text(-0.3 , 1.2, 'C',
#          ha='center',va='center',
#          transform = axB2.transAxes, color = simcolor, fontsize=abcd_font_size)

axB2lower= fig.add_subplot(grid[2,1],xticks=[0,25,50,75,100])
heatmap = axB2lower.pcolor(gamma_arr_Maeda, rho_arr_Maeda, pop_rate_Maeda.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.55)
axB2lower.set_yscale('log')
cbar=fig.colorbar(heatmap, ax=axB2lower);cbar.ax.tick_params(labelsize = tick_font_size) 
axB2lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB2lower.set_title('Maeda 2004\nw/o noise', color= maedacolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

axB2lower.text(-0.3 , 1.2, '(c)',
          ha='center',va='center',
          transform = axB2lower.transAxes, color = simcolor, fontsize=abcd_font_size)

axB3= fig.add_subplot(grid[1,2], xticks=[0,25,50,75,100])
heatmap = axB3.pcolor(k_arr_Gregor_noise, rho_arr_Gregor_noise, pop_rate_Gregor_noise.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,1.2)
axB3.set_yscale('log')
cbar=fig.colorbar(heatmap, ax=axB3);cbar.ax.tick_params(labelsize = tick_font_size) 
axB3.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB3.set_title('Gregor 2010\nwith noise',color= gregorcolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
axB3.text(-0.3 , 1.2, '(d)',
         ha='center',va='center',
         transform = axB3.transAxes, color = simcolor, fontsize=abcd_font_size)

# axB3lower = fig.add_subplot(grid[2,2], xticks=[0,25,50,75,100])
# heatmap = axB3lower.pcolor(k_arr_Gregor, rho_arr_Gregor, pop_rate_Gregor.transpose(), cmap='jet') # cmap='jet'
# heatmap.set_clim(0,1.2)
# axB3lower.set_yscale('log')
# cbar=fig.colorbar(heatmap, ax=axB3lower);cbar.ax.tick_params(labelsize = tick_font_size) 
# axB3lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
# axB3lower.set_title('Gregor 2010\nw/o noise',color= gregorcolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# Sgro regular noise (sig = 0.15)
axB4= fig.add_subplot(grid[1,3],xticks=[0,0.5,1])
heatmap = axB4.pcolor(j_arr_Sgro_regular_noise, rho_arr_Sgro_regular_noise, 
                     pop_rate_Sgro_regular_noise.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.55)
axB4.set_yscale('log'); # axB4.set_ylim([10**(-5),10**(-3)]); 
cbar=fig.colorbar(heatmap, ax=axB4);cbar.ax.tick_params(labelsize = tick_font_size) 
axB4.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB4.set_title('Sgro 2015\nwith noise', color= sgrocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
axB4.text(-0.3 , 1.2, '(e)',
         ha='center',va='center',
         transform = axB4.transAxes, color = simcolor, fontsize=abcd_font_size)
# axB4.set_xticklabels([0,0.25,0.5,0.75,1], rotation=45,fontsize=tick_font_size)

# # Sgro w/o noise  (sig = 0.0)
# axB4lower= fig.add_subplot(grid[2,3],xticks=[0,0.5,1])
# heatmap = axB4lower.pcolor(j_arr_Sgro_no_noise, rho_arr_Sgro_no_noise, 
#                           pop_rate_Sgro_no_noise.transpose(), cmap='jet') # cmap='jet'
# heatmap.set_clim(0,0.55)
# axB4lower.set_yscale('log'); # axB4lower.set_ylim([10**(-5),10**(-3)]); 
# cbar=fig.colorbar(heatmap, ax=axB4lower);cbar.ax.tick_params(labelsize = tick_font_size) 
# axB4lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
# axB4lower.set_title('Sgro 2015\nw/o noise', color= sgrocolor,fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

# # Sgro low noise  (sig = 0.1)
# axB4lower= fig.add_subplot(grid[2,3],xticks=[0,0.25,0.5,0.75,1])
# heatmap = axB4lower.pcolor(j_arr_Sgro_low_noise, rho_arr_Sgro_low_noise, 
#                           pop_rate_Sgro_low_noise.transpose(), cmap='jet') # cmap='jet'
# heatmap.set_clim(0,0.6)
# axB4lower.set_yscale('log'); axB4lower.set_ylim([10**(-5),10**(-3)]); 
# cbar=fig.colorbar(heatmap, ax=axB4lower);cbar.ax.tick_params(labelsize = tick_font_size) 
# axB4lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
# axB4lower.set_title('Sgro 2015\nlow noise', color=mycolors[5],fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})


# axB5= fig.add_subplot(grid[1,4], xticks=[0,25,50,75,100])
# heatmap = axB5.pcolor(gamma_arr_Kamino_noise, rho_arr_Kamino_noise, 
#                      pop_rate_Kamino_noise.transpose(), cmap='jet') # cmap='jet'
# heatmap.set_clim(0,0.65)
# axB5.set_yscale('log')
# cbar=fig.colorbar(heatmap, ax=axB5);cbar.ax.tick_params(labelsize = tick_font_size) 
# axB5.tick_params(axis='both', which='major', labelsize=tick_font_size)
# axB5.set_title('Kamino 2017\nwith noise', color= kaminocolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})
# axB5.text(-0.3 , 1.2, 'F',
#          ha='center',va='center',
#          transform = axB5.transAxes, color = simcolor, fontsize=abcd_font_size)

# # insets
# #ax7= fig.add_axes([0.77,0.13,0.08,0.16])
# axB5in= fig.add_axes([0.835,0.41,0.042,0.082])
# heatmap = axB5in.pcolor(gamma_arr_Kamino_noise, rho_arr_Kamino_noise,
#                        pop_rate_Kamino_noise.transpose(), cmap='jet') # cmap='jet'
# heatmap.set_clim(0,0.65)
# axB5in.set_yscale('log');axB5in.set_xscale('log');
# axB5in.set_xticks([]) ; axB5in.set_yticks([]) 
# axB5in.spines['bottom'].set_color('white');axB5in.spines['top'].set_color('white')
# axB5in.spines['left'].set_color('white');axB5in.spines['right'].set_color('white')

axB5lower= fig.add_subplot(grid[2,4], xticks=[0,25,50,75,100])
heatmap = axB5lower.pcolor(gamma_arr_Kamino, rho_arr_Kamino, 
                          pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.65)
axB5lower.set_yscale('log')
cbar=fig.colorbar(heatmap, ax=axB5lower);cbar.ax.tick_params(labelsize = tick_font_size) 
axB5lower.tick_params(axis='both', which='major', labelsize=tick_font_size)
axB5lower.set_title('Kamino 2017\nw/o noise', color = kaminocolor, fontdict={'fontsize': title_font_size, 'fontweight': 'medium'})

axB5lowerin= fig.add_axes([0.835,0.115,0.042,0.082])
heatmap = axB5lowerin.pcolor(gamma_arr_Kamino, rho_arr_Kamino,
                       pop_rate_Kamino.transpose(), cmap='jet') # cmap='jet'
heatmap.set_clim(0,0.65)
axB5lowerin.set_yscale('log'); axB5lowerin.set_xscale('log');
axB5lowerin.set_xticks([]) ; axB5lowerin.set_yticks([]) 
axB5lowerin.spines['bottom'].set_color('white');axB5lowerin.spines['top'].set_color('white')
axB5lowerin.spines['left'].set_color('white');axB5lowerin.spines['right'].set_color('white')

axB5lower.text(-0.3 , 1.2, '(f)',
          ha='center',va='center',
          transform = axB5lower.transAxes, color = simcolor, fontsize=abcd_font_size)

fig.text(0.51, 0.045, 'Dilution Rate, A.U.',fontsize=label_font_size, ha='center')
fig.text(0.08, 0.35, 'Population Density, A.U.',fontsize=label_font_size, va='center', rotation='vertical')

# plt.subplots_adjust(left=0.1, right=0.96, top=0.93, bottom=0.04, wspace=0.4, hspace=0.72)
# grid.tight_layout(fig5,rect=[0, 0, 1, 1]) 
fig.savefig('fig8_pop_firing_rate_tight_200530.png', bbox_inches='tight') 