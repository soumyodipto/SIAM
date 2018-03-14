import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize







peak_pos_plus=[]
peak_pos_minus=[]
weight_plus=[]
weight_minus=[]






def LinInterpolation(i,j, delE):
    m = (delE[i,j]-delE[i-1,j])/dw
    root=[]
    inc=dw/10.
    for inc in np.arange(dw/10, dw, dw/10):
        root.append(abs(delE[i-1,j]+m*inc))    
    ind = root.index(min(root))
    return omega_grid[i-1]+np.arange(dw/10, dw, dw/10)[ind]

def FindCoeff(i,j, Cn, root):
    m = (Cn[i,j+1]-Cn[i-1,j+1])/dw
    inc = root - omega_grid[i-1]
    return Cn[i-1,j+1]+m*inc






for wmin in np.arange(0., 1., 0.01):
    #print wmin
    Enplus = np.loadtxt('/home/soumyo/Dropbox/Paper1_final/Spectra/nsites300_log_discretization/nimp3/siam_dos_emb_lamda1.05_nsites100_U14_delta1e-08_nimp3_wmin'+str(wmin)+'/En_plus_U')
    Enminus = np.loadtxt('/home/soumyo/Dropbox/Paper1_final/Spectra/nsites300_log_discretization/nimp3/siam_dos_emb_lamda1.05_nsites100_U14_delta1e-08_nimp3_wmin'+str(wmin)+'/En_minus_U')
    Cnplus = np.loadtxt('/home/soumyo/Dropbox/Paper1_final/Spectra/nsites300_log_discretization/nimp3/siam_dos_emb_lamda1.05_nsites100_U14_delta1e-08_nimp3_wmin'+str(wmin)+'/Cn_plus_U')
    Enminus = np.loadtxt('/home/soumyo/Dropbox/Paper1_final/Spectra/nsites300_log_discretization/nimp3/siam_dos_emb_lamda1.05_nsites100_U14_delta1e-08_nimp3_wmin'+str(wmin)+'/Cn_minus_U')
    omega_grid = Enplus[:,0]
    dw = omega_grid[1]-omega_grid[0]
    delEplus = np.zeros([np.shape(Enplus)[0], np.shape(Enplus)[1]-1],float)
    delEminus = np.zeros(np.shape(delEplus),float)
    for i in xrange(np.shape(delEplus)[0]):
        for j in xrange(np.shape(delEplus)[1]):
            delEplus[i,j] = Enplus[i,0] - Enplus[i,j+1]
            delEminus[i,j] = Enminus[i,0] + Enminus[i,j+1]

    for j in xrange(np.shape(delEplus)[1]):
        for i in xrange(len(omega_grid)):
            if i > 0:
                if delEplus[i,j]*delEplus[i-1,j] < 0:
                    root = LinInterpolation(i,j,delEplus)
                    peak_pos_plus.append(root)
                    grad_plus = (Enplus[i,j+1] - Enplus[i-1,j+1])/dw
                    coeff_plus = FindCoeff(i, j, Cnplus, root)
                    weight_plus.append(coeff_plus/(abs(1.-grad_plus)))
              
	        if delEminus[i,j]*delEminus[i-1,j] < 0:
		    root = LinInterpolation(i,j,delEminus)
		    peak_pos_minus.append(root)	
		    grad_minus = (Enminus[i,j+1] - Enminus[i-1,j+1])/dw
		    coeff_minus = FindCoeff(i, j, Cnminus, root)
		    weight_minus.append(coeff_minus/(abs(1.+grad_minus)))








omega_min = min(peak_pos_plus)


for omega_min in peak_pos_plus:
    if omega_min > 0:
        break



weight = weight_minus + weight_plus
weight += weight
print len(weight)
#weight = 1/sum(weight)*np.array(weight)


#for ind,height in enumerate(weight):
#    if height >0.07:
#        weight[ind]=0.

#print sum(weight)
peak_pos = peak_pos_minus + peak_pos_plus
peak_pos += list(-1.*np.array(peak_pos))
#weight = 1/sum(weight)*np.array(weight)
print sum(weight)

#plt.plot(omega_grid, Enplus[:,1]-omega_grid, 'b-')
#plt.plot(omega_grid, Cnplus[:,1], 'b-')
#plt.axhline(0, color='black')
#plt.show()


plt.plot(peak_pos, weight, 'ro')
plt.xlim(-1,1)
#plt.ylim(0,0.3)
plt.vlines(peak_pos, [0], weight)
plt.show()


#fd=file("siam_nonint_dos_nsites"+str(300)+'_nimp'+str(1)+'_lambda'+str(2.0),"w")
#for i, height in enumerate(weight):
#   print >>fd, peak_pos[i], weight[i]


'''
mom1 = 0.
for i in xrange(len(peak_pos)):
   mom1 += weight[i]*(peak_pos[i]**2)

print mom1
'''

def LorentzianBroaden(peak_pos):
    eta = 0.05
    spectral_fn = []
    for omega in np.arange(-1,1, 0.001):
        gf = 0.
        for i, peak in enumerate(peak_pos):
            gf += 1/np.pi*eta*weight[i]*1/((omega-peak)**2+(eta)**2)
        spectral_fn.append(gf)
    return spectral_fn


def GaussBroaden(omega, peak, c1, c2):
    b = c1*omega_min + c2*abs(peak)
    #b = omega_min + 1.5*abs(peak)
    #b = 0.05
    return 1./(b*np.sqrt(np.pi))*np.exp(-(((omega-peak)/b)**2))

def MakeGaussBroadenedSpectra(peak_pos, c1, c2):
    spectral_fn=[]    
    for omega in np.arange(-1,1, 0.001):
        gf = 0.
        for i, peak in enumerate(peak_pos):
            gf += weight[i]*GaussBroaden(omega, peak, c1, c2)
        spectral_fn.append(gf)
    return spectral_fn



def FindHubBandMax(spectral_fn):
    for i in xrange(1, len(spectral_fn)-1):
        if spectral_fn[i] > spectral_fn[i-1] and spectral_fn[i] > spectral_fn[i+1]:
            hub_max = spectral_fn[i]
            hub_max_en = np.arange(-.1,.1,0.001)[i]
            break
    return hub_max_en, hub_max
    


def FindDosFermiEnergy(peak_pos, weight, c1, c2):
    s = 0.
    for i in xrange(len(peak_pos)):
        s += 1./((c1*omega_min+c2*abs(peak_pos[i]))*np.sqrt(np.pi))*weight[i]*np.exp(-(peak_pos[i]**2)/(c1*omega_min+c2*abs(peak_pos[i]))**2)
    return s



'''
c2=0.35
for c1 in np.arange(10.5,11.5,0.05):
    spectral_fn = MakeGaussBroadenedSpectra(peak_pos, c1, c2)    
    area = 0.
    for num in xrange(len(np.arange(-1.,1.,0.001))-1):
        area += ((spectral_fn[num]+spectral_fn[num+1])*abs(dw))
    area*=0.5
    print c1, 0.05*np.pi/area*np.array(spectral_fn)[1000]
'''





spectral_fn = MakeGaussBroadenedSpectra(peak_pos, 1., 0.4)
#spectral_fn = LorentzianBroaden(peak_pos)

area = 0.
for num in xrange(len(np.arange(-1.,1.,0.001))-1):
    area += ((spectral_fn[num]+spectral_fn[num+1])*abs(dw))

area*=0.5    
print area


plt.plot(np.arange(-1.,1.,0.001), 0.05*np.pi/area*np.array(spectral_fn), 'b-')

#nonintdos = np.loadtxt('siam_dos_nsites300_exact_U0')

#plt.plot(np.arange(-1.,1.,0.001), 0.05*np.pi*nonintdos[:,1], 'r-')

#plt.legend(('Emb(2)', 'Exact'))
plt.xlim(-1,1)
plt.ylim(0,1)
plt.show()




'''
fd = open('siam_dos_nsites300_nimp1_U2','w')
for i in xrange(len(spectral_fn)):
    print >>fd, np.arange(-1.,1.,0.001)[i],  0.05*np.pi/area*np.array(spectral_fn)[i]
    
'''


'''
fd = open('siam_dos_peak_nsites300_nimp2_U10_delta1e-12','w')
for i in xrange(len(peak_pos)):
    print >>fd, peak_pos[i],  weight[i]
'''
