from os import path
import h5py
import numpy as N
import matplotlib
import matplotlib.patches as mpatches
import pylab as plt


plt.rc('text.latex', preamble = '\usepackage{amsmath},' '\usepackage{yfonts},' '\usepackage[T1]{fontenc},' '\usepackage[latin1]{inputenc},' '\usepackage{txfonts},' '\usepackage{times},' '\usepackage{blindtext}' )
       

       
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':8})
plt.rc('lines', linewidth=1.5)
plt.rc('xtick', labelsize='medium')
plt.rc('ytick', labelsize='medium')
plt.rc('legend', fontsize='medium')
plt.rc('text', usetex=True)


exact=N.loadtxt("siam_ed_gs_nsites8",usecols=(0,1))
uhf=N.loadtxt("siam_uhf_gs_nsites8",usecols=(0,1))
rhf=N.loadtxt("siam_rhf_gs_nsites8",usecols=(0,1))
emb1=N.loadtxt("siam_emb_gs_nsites8_nimp1",usecols=(0,1))
emb2=N.loadtxt("siam_emb_gs_nsites8_nimp2",usecols=(0,1))
emb3=N.loadtxt("siam_emb_gs_nsites8_nimp3",usecols=(0,1))
emb4=N.loadtxt("siam_emb_gs_nsites8_nimp4",usecols=(0,1))
#var_emb=N.loadtxt("VariationalEmbedding-Siam-8-sites",usecols=(0,1))

exactx=exact[:,0]/0.1
exacty=exact[:,1]

uhfx=uhf[:,0]/0.1
uhfy=uhf[:,1]

rhfx=rhf[:,0]/0.1
rhfy=rhf[:,1]


embx=emb1[:,0]/0.1
emby=emb1[:,1]

emb2x=emb2[:,0]/0.1
emb2y=emb2[:,1]

emb3x=emb3[:,0]/0.1
emb3y=emb3[:,1]

emb4x=emb4[:,0]/0.1
emb4y=emb4[:,1]


fig_size = (3.375,3.375/1.61)

f1, ax1 = plt.subplots(1, 1, figsize=fig_size, dpi=200)

#col = 'black'

line0, =ax1.plot(exactx,exacty,'k')
line1, =ax1.plot(uhfx,uhfy,'g',markevery=2)
line2, =ax1.plot(rhfx,rhfy,'b',markevery=2)
line3, =ax1.plot(embx,emby,'c',markevery=2)
line4, =ax1.plot(emb2x,emb2y,'m--')
line5, =ax1.plot(emb4x,emb4y,'ro', markevery=4, markersize=3.5)

col = 'black'
ax1.legend(('ED','UHF','RHF','Emb(1)', 'Emb(2)','Emb(4)'),ncol=1,loc='upper right', frameon=True, prop={'size':6.3})
ax1.xaxis.set_ticks(N.arange(0,12,step=2))
ax1.yaxis.set_ticks(N.arange(-0.59,-0.51,step=0.02))
#ax1.yaxis.set_ticks([0.0,0.5,1.0,1.5,2.0])
ax1.set_xlabel('$U/V$')
ax1.set_ylabel('$E_{gs}/L$')


ax1.set_xlim(0,10)
ax1.set_ylim(-0.59, -0.53)
plt.grid(True)
f1.tight_layout(pad=0.15)
f1.subplots_adjust(wspace=0)


f1.savefig('fig_1.pdf')
plt.show()
