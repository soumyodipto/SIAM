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
import scipy.optimize
import embed
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
    #N=[]
    #N_s=[]
    #N_b=[]
    for occ_list in product(occ, repeat=L):
        state=MakeFock(occ_list)        
        state_dw = state[:L]
        state_up = state[L:]
        Ndw = sum_digits(int(state_dw))
	Nup = sum_digits(int(state_up))
	if Nup+Ndw==L and Nup==Ndw:
	#if Nup==2 and Ndw==3:
	    fock.append(state)
	   
        #N.append(sum_digits(int(state)))        
        #fock.append(state)
    fock.remove('')
    return fock


def AddCoreVirtToFock(nimp, N_tot, Nup, Ndw, fock):
    ncore_up = Nup-nimp
    nvirt_up = N_tot - 2*nimp - ncore_up
    ncore_dw = Ndw-nimp
    nvirt_dw = N_tot - 2*nimp - ncore_dw
    corevirt_up = ['1']*ncore_up + ['0']*nvirt_up
    #corevirt_up = "".join(corevirt_up)
    corevirt_dw = ['1']*ncore_dw + ['0']*nvirt_dw
    #corevirt_dw = "".join(corevirt_dw)
    fock_full=[]
    for i in fock:
        fock_up=list(i[:2*nimp])
        fock_dw=list(i[2*nimp:])
        fock_up=fock_up+corevirt_up
        fock_dw=fock_dw+corevirt_dw
        full=fock_up+fock_dw
        full="".join(full)
        fock_full.append(full)
    return fock_full



def Swap(basis,j,k):
    new = basis[:j]+basis[k]+basis[j+1:k]+basis[j]+basis[k+1:]
    if k-j==1:
        phase=1.
    else:    
        phase = (-1.)**(sum_digits(int(basis[j+1:k])))
    return new, phase    
    
    


def MakeCorrelatedHam(nimp,l,L,fock,U,hup,hdw):
    H_U=utils.zeros([l,l])
    H_diag=utils.zeros([l,l])
    H_t=utils.zeros([l,l])
    for i in xrange(len(fock)):
        basis=fock[i]
        for j in xrange(L):
            if basis[j]==basis[j+L] and basis[j]=='1':
                if j == 0:
                    H_U[i,i]+=U
                H_diag[i,i]+=hup[j,j]
                H_diag[i,i]+=hdw[j,j]
            if basis[j] != basis[j+L]:
                if basis[j] == '1':
                    H_diag[i,i]+=hup[j,j]
                else:
                    H_diag[i,i]+=hdw[j,j]                          
	    for k in xrange(L):
		if j < k:
		    if basis[j] != basis[k]:
		        new,phase=Swap(basis,j,k)
		        if new in fock:
		            ind = fock.index(new)
		            H_t[i,ind]=phase*hup[j,k]
		    if basis[j+L] != basis[k+L]:
		        new,phase=Swap(basis,j+L,k+L)
		        if new in fock:
		            ind = fock.index(new)
		            H_t[i,ind]=phase*hdw[j,k]    
    return H_U,H_diag,H_t    
    
    
    
    
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
   


def SIAMEmbGS(nimp, nsites, U, V, cb_edge):

    U *= V
    N_tot=nsites
    Nup=N_tot/2
    Ndw=N_tot-Nup
    
    hlat = MakeSingleParticleHam(nsites, V, cb_edge) #defined the single particle hopping energies
    fock = GenerateFockSpace(2*nimp)
    fockfull = AddCoreVirtToFock(nimp, N_tot, Nup, Ndw, fock)
    l=len(fock)


    hlat[0,0]=-0.5*U
    fock_mat = 1.*hlat
    fock_mat[0,0] += 0.5*U
    pemb, pcore, pvirt = embed.embed(fock_mat, nimp, Nup)
    basis = N.hstack([pemb, pcore, pvirt])
    h1e_emb = mdot(N.transpose(basis),hlat,basis)
    H_corr_U, H_corr_diag, H_corr_hop = MakeCorrelatedHam(nimp,l,nsites,fockfull,U,h1e_emb,h1e_emb)
    H_corr = H_corr_U + H_corr_diag + H_corr_hop

    if nimp == 1:
        e,v = numpy.linalg.eigh(H_corr)
        psi_corr = v[:,0]
             
    else:
        e,v=scipy.sparse.linalg.eigsh(scipy.sparse.csr_matrix(H_corr),k=1,which='SA',maxiter=10000,tol=1e-6)
	psi_corr = v[:,0]

    return e[0], mdot(N.transpose(psi_corr), H_corr_U, psi_corr)/U
    

    
def main(nimp, nsites, U, V, cb_edge):

    en, docc = SIAMEmbGS(nimp, nsites, U, V, cb_edge) #calculates the GS energy and avg. double occupation
    print "Number of sites (including the impurity) = %d \n" %nsites

    print "Interaction U = %s \n" %U
    print "Impurity to conduction band hopping V = %s \n" %V
    print "Conduction band edge = %s \n" %cb_edge
    
    print "SIAM ground state energy = %s \n" %en
    print "SIAM average double occupation = %s \n" %docc

    return 0.
    

    

if __name__ == "__main__":
    main(1, 8, 2., 1., 1.)    
