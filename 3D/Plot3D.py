#Plotting data for 3D Ising Model 
#Columns from left: separation, avg random, sd random, avg order, sd order, 
#avg hilbert, sd hilbert, avg lebesgue, sd lebesgue

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from matplotlib import style

style.use( 'ggplot' ) #style of plot

beta = 0.6 #temperature
N = 16*16*16 #total number of sites (3d square lattice of side 16)
MCS = 10000 #total number of updates
BS = 100 #bin size

#load the data from a file Data_3D_*temp with 2 decimal places*.csv
df = pd.read_csv( 'Data_3D_%.2f.csv' % beta, skiprows=0, header=1 )

#cols 1*n are for avg.s and cols 2*n are for s.d.s
mean = np.zeros( (11, 8) ) 
std = np.zeros( (11, 8) ) 

for i in range( 0, 11, 1 ): #looping over seperation
    temp_mean = df.loc[df['separation'] == i].mean( axis=0 )[1:]
    temp_std = df.loc[df['separation'] == i].std( axis=0 )[1:]
    
    for j in range( 0, 8, 1 ): #looping over methods
        mean[i][j] = temp_mean[j]
        std[i][j] = temp_std[j]

x = np.arange( 0, 11, 1 )

labels = [ 'Random', 'Order', 'Hilbert', 'Lebesgue' ]

fig, ax = plt.subplots( nrows=1, ncols=2, sharex=True, figsize=(17, 6) )

for i in range( 0, 4, 1 ):
    ax[0].errorbar( x, mean[:, 2*i], yerr=std[:, 2*i], capsize=3, 
                   label=labels[i] )
    ax[1].errorbar( x, mean[:, 2*i+1], yerr=std[:, 2*i+1], capsize=3, 
                   label=labels[i] )

ax[0].set_xlabel( 'Separation' )
ax[0].set_ylabel( 'Mean Correlation' )
ax[0].set_title( '3D Ising Model \nCorrelation vs. Separation \n'
                +r'$(\beta=%.2f, N=%d, MCS=%d, BS=%d)$' % (beta, N, MCS, BS) )

ax[1].set_xlabel( 'Separation' )
ax[1].set_ylabel( 'Standard Deviation' )
ax[1].set_title( '3D Ising Model \nStandard Deviation of  Mean Correlation '
                +'vs. Separation \n'
                +r'$(\beta=%.2f, N=%d, MCS=%d, BS=%d)$' % (beta, N, MCS, BS) )

ax[0].legend()
ax[1].legend()

plt.tight_layout()
plt.show()





















