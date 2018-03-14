from os import path
import h5py
import numpy as np
import matplotlib
import matplotlib.patches as mpatches
import pylab as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import FormatStrFormatter


majorFormatter = FormatStrFormatter('%g')

plt.rc('text.latex', preamble = '\usepackage{amsmath},' '\usepackage{yfonts},' '\usepackage[T1]{fontenc},' '\usepackage[latin1]{inputenc},' '\usepackage{txfonts},' '\usepackage{times},' '\usepackage{blindtext},' '\usepackage{braket}' )
       

       
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':8})
plt.rc('lines', linewidth=0.7)
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
plt.rc('legend', fontsize='medium')
plt.rc('text', usetex=True)

nsites=8
delta=0.005
nimps=[1,2]
Us = [0.0, 2.5, 5.0, 7.5]

omega_grid = np.arange(-1.,1.,0.01)


edy=[]
emb1y=[]
emb2y=[]
#emb4y=[]
for ind, U in enumerate(Us):
    eddir = './siam_ED/siam_dos_ed_U'+str(U)+'_nsites'+str(nsites)+'_delta'+str(delta)+'/siam_dos_ed_nsites'+str(nsites)+'_U'+str(U)+'_eta'+str(delta)
    data_ed = np.loadtxt(eddir)
    edy.append(data_ed[:,1])
    emb1dir = './siam_emb/siam_emb_dos_nsites'+str(nsites)+'_U'+str(U/10)+'_nimp'+str(1)+'_eta'+str(delta)
    data_emb1 = np.loadtxt(emb1dir)
    emb1y.append(data_emb1[:,1])
    emb2dir = './siam_emb/siam_emb_dos_nsites'+str(nsites)+'_U'+str(U/10)+'_nimp'+str(2)+'_eta'+str(delta)
    data_emb2 = np.loadtxt(emb2dir)
    emb2y.append(data_emb2[:,1])
    #emb4dir = './siam_emb/siam_emb_dos_nsites'+str(nsites)+'_U'+str(U/10)+'_nimp'+str(4)+'_eta'+str(delta)
    #data_emb4 = np.loadtxt(emb4dir)
    #emb4y.append(data_emb4[:,1])


fig_size = (3.375,3.375/1.3)
f1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=fig_size, dpi=200, sharex='col', sharey='row')

#dashes = [5,2,10,5]
dashes=[2, 1]
dashes1 = [5,2] 
#ax1 = plt.subplot(221)
line1, = ax1.plot(omega_grid, edy[0],'k', markersize=1.5,markevery=1)
line2, = ax1.plot(omega_grid, emb1y[0], 'b--')
line3, = ax1.plot(omega_grid, emb2y[0], 'r--')
line2.set_dashes(dashes)
line3.set_dashes(dashes1)
#line4, = ax1.plot(omega_grid, emb4y[0], 'r--')
#ax1.legend(('ED','Emb(2)','Emb(1)'),'upper right',ncol=1,prop={'size':6})
ax1.set_ylim(0,25)
#textstr = '$\eta='+str(delta)+'$'
#plt.text(0.1, 0.7, textstr ,verticalalignment='top', horizontalalignment='left',transform=ax1.transAxes, color='black', fontsize=20)
textstr = '$U/V=0$'
plt.text(0.1, 0.9, textstr ,verticalalignment='top', horizontalalignment='left',transform=ax1.transAxes, color='black')
#ax1.set_xlabel(r'$\omega$')
ax1.set_ylabel(r'$A(\omega)$')
ax1.yaxis.set_ticks(np.arange(0,25,10))
minorLocator = AutoMinorLocator()
ax1.yaxis.set_minor_locator(minorLocator)
ax1.tick_params(which='both', width=0.5)
ax1.tick_params(which='major', length=4)
ax1.tick_params(which='minor', length=1.5)
ax1.yaxis.set_major_formatter(majorFormatter)




#ax2 = plt.subplot(222, sharex=ax1, sharey=ax1)
line1, = ax2.plot(omega_grid, edy[1],'k', markersize=1.5,markevery=1)
line2, = ax2.plot(omega_grid, emb1y[1], 'b--')
line3, = ax2.plot(omega_grid, emb2y[1], 'r--')
line2.set_dashes(dashes)
line3.set_dashes(dashes1)
ax2.legend(('ED','Emb(1)','Emb(2)'),'upper right',frameon=True,ncol=1,prop={'size':6})
#line4, = ax2.plot(omega_grid, emb4y[1], 'r--')
textstr = '$U/V=2.5$'
plt.text(0.1, 0.9, textstr ,verticalalignment='top', horizontalalignment='left',transform=ax2.transAxes, color='black')
#ax2.tick_params(axis='x',length=0, width=0)



#ax3 = plt.subplot(223, sharex=ax1, sharey=ax1)
line1, = ax3.plot(omega_grid, edy[2],'k', markersize=1.5,markevery=1)
line2, = ax3.plot(omega_grid, emb1y[2], 'b--')
line3, = ax3.plot(omega_grid, emb2y[2], 'r--')
line2.set_dashes(dashes)
line3.set_dashes(dashes1)
#line4, = ax3.plot(omega_grid, emb4y[2], 'r--')
ax3.set_xlabel(r'$\omega$')
ax3.set_ylabel(r'$A(\omega)$')
ax3.set_ylim(0,25)
textstr = '$U/V=5.0$'
plt.text(0.1, 0.9, textstr ,verticalalignment='top', horizontalalignment='left',transform=ax3.transAxes, color='black')
ax3.yaxis.set_ticks(np.arange(0,25,10))
ax3.xaxis.set_ticks([-1,-0.5,0,0.5,1])
minorLocator = AutoMinorLocator()
ax3.yaxis.set_minor_locator(minorLocator)
ax3.tick_params(which='both', width=0.5)
ax3.tick_params(which='major', length=4)
ax3.tick_params(which='minor', length=1.5)
ax3.xaxis.set_major_formatter(majorFormatter)
ax3.yaxis.set_major_formatter(majorFormatter)


#ax4 = plt.subplot(224, sharex=ax1, sharey=ax1)
line1, = ax4.plot(omega_grid, edy[3],'k', markersize=1.5,markevery=1)
line2, = ax4.plot(omega_grid, emb1y[3], 'b--')
line3, = ax4.plot(omega_grid, emb2y[3], 'r--')
line2.set_dashes(dashes)
line3.set_dashes(dashes1)
#line4, = ax4.plot(omega_grid, emb4y[3], 'r--')
textstr = '$U/V=7.5$'
plt.text(0.1, 0.9, textstr ,verticalalignment='top', horizontalalignment='left',transform=ax4.transAxes, color='black')
ax4.set_xlabel(r'$\omega$')
ax4.xaxis.set_major_formatter(majorFormatter)



f1.subplots_adjust(hspace=0)
#f.subplots_adjust(wspace=0.05)
f1.tight_layout(pad=0.15)
f1.subplots_adjust(wspace=0.1)


f1.savefig('fig_3.pdf')
plt.show()

