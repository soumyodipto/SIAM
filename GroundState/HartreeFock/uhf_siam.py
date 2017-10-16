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


def ReturnSiamUHFFockMat(N_tot, Nup, Ndw, hlat, U):
	nup_avg=0.
	ndw_avg=1.
	uhf_convg=False
	for it in xrange(500):
	    hlat_up = 1.*hlat
	    hlat_up[0,0] += U*(ndw_avg)
	    hlat_dw = 1.*hlat
	    hlat_dw[0,0] += U*(nup_avg)
	    e_up,v_up = numpy.linalg.eigh(hlat_up)
	    e_dw,v_dw = numpy.linalg.eigh(hlat_dw)
	    nup=0.
	    ndw=0.
	    for k in xrange(Nup):
		nup+=(v_up[0,k]**2)
	    for k in xrange(Ndw):
		ndw+=(v_dw[0,k]**2)    
	    if abs(nup-nup_avg) < 1e-8 and abs(ndw-ndw_avg) < 1e-8:
		print "\nUHF SIAM Converged for U = %s" % U
		en=0.
		for k in xrange(Nup):
		    en+=e_up[k]
		for k in xrange(Ndw):
		    en+=e_dw[k]
		en-=(U*nup*ndw)
		uhf_convg=True    
		break
	    else:
		nup_avg=nup
		ndw_avg=ndw

	if uhf_convg:
	    return uhf_convg, en, hlat_up, hlat_dw
	else:
	    print "UHF Siam did not converge for U = %s" % U
            return uhf_convg, None, None, None
            


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




# Main function that solves for the SIAM ground state using Unestricted Hartree Fock
# SIAM input parameters are defined - nsites, U, V, cbedge
def main(nsites, V, cb_edge):
    hlat = MakeSingleParticleHam(nsites, V, cb_edge) #defined the single particle hopping energies
    hlat[0,0]=-0.5*U
    N_tot=nsites #Number of electrons (half-filling in this case)
    Nup=N_tot/2 #Number of up spin electrons
    Ndw=N_tot-Nup #Number of down spin electrons
    uhf_convg, en, hlat_up, hlat_dw = ReturnSiamUHFFockMat(N_tot, Nup, Ndw, hlat, U)
    
    print "Number of sites (including the impurity) = %d \n" %nsites

    print "Interaction U = %s \n" %U
    print "Impurity to conduction band hopping V = %s \n" %V
    print "Conduction band edge = %s \n" %cb_edge
    
    #print "SIAM ground state energy (Restricted Hartree Fock) = %s \n" %en

    #SIAM ground state energy vs U(Unrestricted Hartree Fock)
    fd=file('siam_uhf_gs_nsites'+str(nsites),"w")

    #SIAM ground state double occupation vs U(Unrestricted Hartree Fock)
    fd1=file('siam_uhf_gs_docc_nsites'+str(nsites),"w")
    
    for U in N.arange(0.1,10.1,0.1):  #Scan of the interaction U
        U*=V
        hlat[0,0]=-0.5*U
        uhf_convg, en, hlat_up, hlat_dw = ReturnSiamUHFFockMat(N_tot, Nup, Ndw, hlat, U)
        #Find the double occupation of the uhf fock matrices 
        if uhf_convg:
	    print >>fd, U, en/N_tot
	    e_up,v_up = numpy.linalg.eigh(hlat_up)
	    e_dw,v_dw = numpy.linalg.eigh(hlat_dw)
	    nup=0.
	    ndw=0.
	    for k in xrange(Nup):
		nup+=(v_up[0,k]**2)
	    for k in xrange(Ndw):
		ndw+=(v_dw[0,k]**2)
	    print >>fd1, U, nup*ndw

    return 0
    
#Calls the main function
if __name__ == "__main__":
    main(6, 1., 1.)     
