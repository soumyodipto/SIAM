#!/usr/bin/env python
import math
import time
import numpy as N
import numpy.linalg
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import scipy.optimize


#Returns rhf convergence, rhf energy
def ReturnSiamRHFFockMat(N_tot, Nup, Ndw, hlat, U):
	nup_avg=0.5
	rhf_convg=False
	for it in xrange(500):
	    hlat_up = 1.*hlat
	    hlat_up[0,0] += U*(nup_avg)
	    hlat_dw = 1.*hlat
	    hlat_dw[0,0] += U*(nup_avg)
	    e_up,v_up = numpy.linalg.eigh(hlat_up)
	    e_dw,v_dw = numpy.linalg.eigh(hlat_dw)
	    nup=0.
	    for k in xrange(Nup):
		nup+=(v_up[0,k]**2)    
	    if abs(nup-nup_avg) < 1e-8 :
		print "\nRHF SIAM Converged for U = %s" % U
		en=0.
		for k in xrange(Nup):
		    en+=e_up[k]
		en*=2.
		en-=(U*nup*nup)
		rhf_convg=True    
		break
	    else:
		nup_avg=nup

	if rhf_convg:
	    return rhf_convg, en
	else:
	    print "RHF Siam did not converge for U = %s" % U
            return rhf_convg, None


#Defines the single particle hamiltonian or the single particle hopping energies
def MakeSingleParticleHam(nsites, V, cb_edge):
    hlat=utils.zeros([nsites,nsites]) #the indexing of hlat is imp, cbsite1, -cbsite1, ... , cbsite(nsites-2)/2,-cbsite(nsites-2)/2 
    hlat[0,1:] = V #ASSUMES a constant impurity to conduction band hopiing
    hlat[1:,0] = V
    cb_states = [0] #list of conduction band energy levels
    for i in range((nsites-2)/2):
	cb_states.append((i+1)/((nsites-2)/2))
	cb_states.append(-(i+1)/((nsites-2)/2))
    for i in xrange(nsites-1):
        hlat[i+1,i+1] = cb_states[i]
    return hlat




# Main function that solves for the SIAM ground state
# SIAM input parameters are defined - nsites, U, V, cbedge
def main(nsites, U, V, cb_edge):
    hlat = MakeSingleParticleHam(nsites, V, cb_edge) #defined the single particle hopping energies
    hlat[0,0]=-0.5*U
    N_tot=nsites #Number of electrons (half-filling in this case)
    Nup=N_tot/2 #Number of up spin electrons
    Ndw=N_tot-Nup #Number of down spin electrons
    rhf_convg, en = ReturnSiamRHFFockMat(N_tot, Nup, Ndw, hlat, U)
    
    print "Number of sites (including the impurity) = %d \n" %nsites

    print "Interaction U = %s \n" %U
    print "Impurity to conduction band hopping V = %s \n" %V
    print "Conduction band edge = %s \n" %cb_edge
    
    print "SIAM ground state energy (Restricted Hartree Fock) = %s \n" %en
    return 0
    
#Calls the main function
if __name__ == "__main__":
    main(6, 4., 1., 1.)
