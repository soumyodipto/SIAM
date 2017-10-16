#!/usr/bin/env python
import math
import time
import numpy as N
import numpy.linalg
import scipy.linalg
import utils
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from itertools import *




def sum_digits(n):
    s = 0
    while n:
        s += n % 10
        n /= 10
    return s


def mdot(*args):
    """mdot(A,B,C,..) form matrix product A*B*C*..."""
    r = args[0]
    for a in args[1:]:
        r = N.dot(r,a)
    return r






def MakeFock(occ_list):
    num_sites=len(occ_list)
    fockup='0'*(num_sites)
    fockdw='0'*(num_sites)
    fockup=list(fockup)
    fockdw=list(fockdw)
    for i in xrange(num_sites):
        if occ_list[i] == 'down':
            fockdw[i]='1'
        if occ_list[i] == 'up':
            fockup[i] = '1'
        if occ_list[i] == 'updown':
            fockup[i] = '1'
            fockdw[i]='1'
    fock="".join(fockdw+fockup)
    return fock


def GenerateFockSpace(L):
    occ=['vac','down','up','updown']
    fock=['']
    for occ_list in product(occ, repeat=L):
        state=MakeFock(occ_list)        
        state_dw = state[:L]
        state_up = state[L:]
        Ndw = sum_digits(int(state_dw))
	Nup = sum_digits(int(state_up))
	if Nup+Ndw==L and Nup==Ndw:
	    fock.append(state)
    fock.remove('')
    return fock






def Swap(basis,j,k):
    new = basis[:j]+basis[k]+basis[j+1:k]+basis[j]+basis[k+1:]
    if k-j==1:
        phase=1.
    else:    
        phase = (-1.)**(sum_digits(int(basis[j+1:k])))
    return new, phase   



def MakeHamManyBodyBasis(l,L,fock,U,h):
    H_U=utils.zeros([l,l]) #interaction
    H_diag=utils.zeros([l,l]) #on-site energies
    H_t=utils.zeros([l,l]) #inter-site hopping
    for i in xrange(len(fock)):
        basis=fock[i]
        for j in xrange(L):
            if basis[j]==basis[j+L] and basis[j]=='1':
                if j == 0:
                    H_U[i,i]+=U
                H_diag[i,i]+=2.*h[j,j]
            if basis[j] != basis[j+L]:
                H_diag[i,i]+=h[j,j]                          
	    for k in xrange(L):
		if j < k:
		    if basis[j] != basis[k]:
		        new,phase=Swap(basis,j,k)
		        ind = fock.index(new)
		        H_t[i,ind]=phase*h[j,k]
		    if basis[j+L] != basis[k+L]:
		        new,phase=Swap(basis,j+L,k+L)
		        ind = fock.index(new)
		        H_t[i,ind]=phase*h[j,k]    
    return H_U,H_diag,H_t
    
    
    
def SolveExactGSSiam(nsites, hlat, U):
    fock = GenerateFockSpace(nsites) #makes the half-filling basis for nsites in second quantization
    l=len(fock) #size of the basis
    H_U, H_diag, H_hop = MakeHamManyBodyBasis(l,nsites,fock,U,hlat) #defines the ham in the full basis - the interaction part, diagonal part from on site energies, and the ham from hopping 
    H = H_U+H_diag+H_hop #full Hamiltonian
    #Diagonalize using sparse linear algebra routine 
    e,v=scipy.sparse.linalg.eigsh(scipy.sparse.csr_matrix(H),k=1,which='SA',maxiter=10000,tol=1e-6)
    psi0=v[:,0] #ground state eigenvector
    en=e[0] #ground state energy
    docc=mdot(N.transpose(psi0),H_U,psi0)/U #avg. double occupation 
    return en, docc


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
    en, docc = SolveExactGSSiam(nsites, hlat, U) #calculates the GS energy and avg. double occupation

    print "Number of sites (including the impurity) = %d \n" %nsites

    print "Interaction U = %s \n" %U
    print "Impurity to conduction band hopping V = %s \n" %V
    print "Conduction band edge = %s \n" %cb_edge
    
    print "SIAM ground state energy = %s \n" %en
    print "SIAM average double occupation = %s \n" %docc

    return 0.


#Calls the main function
if __name__ == "__main__":
    main(6, 4., 1., 1.)




    
    
