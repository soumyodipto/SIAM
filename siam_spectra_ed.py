#!/usr/bin/env python
import math
import time
import numpy as N
import numpy.linalg
import scipy.linalg
import utils
import scipy.sparse
import scipy.sparse.linalg





def mdot(*args):
    """mdot(A,B,C,..) form matrix product A*B*C*..."""
    r = args[0]
    for a in args[1:]:
        r = N.dot(r,a)
    return r


def sum_digits(n):
    s = 0
    while n:
        s += n % 10
        n /= 10
    return s


def GenerateFockSpace(L,plus,minus):
    fock = ['']
    for J in xrange(4**L):
    	a = bin(J)[2:]
    	b = 2*L - len(a)
    	state = '0'*b + a
    	state_dw = state[:L]
    	state_up = state[L:]
    	Ndw = sum_digits(int(state_dw))
    	Nup = sum_digits(int(state_up))
	if plus:
            if Nup+Ndw==L+1 and abs(Nup-Ndw)==1:
		fock.append(state)
	elif minus:
	    if Nup+Ndw==L-1 and abs(Nup-Ndw)==1:
		fock.append(state)
	else:    	
	    if Nup+Ndw==L and Nup==Ndw:
    	        fock.append(state)    
    fock.remove('')
    return fock


def Swap(basis,j,k):
    new = basis[:j]+basis[k]+basis[j+1:k]+basis[j]+basis[k+1:]
    return new
    

def MakeHerm(mat):
    mat_herm=utils.zeros(mat.shape)
    mat_conj=mat.conjugate()  
    for i in xrange(mat.shape[0]): 
        for j in xrange(mat.shape[0]):
            if i == j:
		mat_herm[i,i]=mat[i,i]
	    elif j > i:
                mat_herm[i,j]=mat[i,j]
	    else:
                mat_herm[i,j]=mat_conj[i,j]
    return mat_herm
    
def cbdiagelements(cb_edge,L,N):
    if L == 2:
        return cb_edge
    N_cb = L - 1
    gap = 2.0*cb_edge/(N_cb-1)
    return -cb_edge + N*gap

def MakeHam(l,L,fock,U,V,cb_edge):
	H_imp = utils.zeros([l,l])
	H_cb = utils.zeros([l,l])
	H_hyb = utils.zeros([l,l])
	#H_U = zeros((l,l),float)
	for i in xrange(l):
	    basis = fock[i]
	    if basis[0] != basis[L]:
		H_imp[i,i] = -U/2
	    #if basis[0] == basis[L] and basis[0] == '1':
		#H_U[i,i] = 1.
	    for j in range(1,L):
		if basis[j] == basis[j+L] and basis[j] == '1':
		    H_cb[i,i] += 2.*cbdiagelements(cb_edge,L,j-1)
		if basis[j] != basis[j+L]:
		    H_cb[i,i] += cbdiagelements(cb_edge,L,j-1)

		if basis[0] != basis[j]:
		    new = basis[j] + basis[1:j] + basis[0] + basis[(j+1):]
		    if j == 1:
			phase = 1.
		    else:
			phase = (-1.)**(sum_digits(int(basis[1:j])))
		    if new in fock:
		        ind = fock.index(new)
		        H_hyb[i,ind] = phase*V
		if basis[L] != basis[j+L]:
		    new = basis[:L] + basis[j+L] + basis[L+1:j+L] + basis[L] + basis[(j+L+1):]
		    if j == 1:
			phase = 1.
		    else:
			phase = (-1.)**(sum_digits(int(basis[L+1:j+L])))
		    if new in fock:	
		        ind = fock.index(new)
		        H_hyb[i,ind] = phase*V

	H = H_imp + H_cb + H_hyb
	return H


def Evaluate_VPsi0(psi0,fock,fock_minus,fock_plus,plus,minus):
    if minus:        
        coeff=utils.zeros(len(fock_minus))        
        for i in xrange(len(fock)):
            bas=fock[i]
            if bas[0] == '1':
                new = '0'+bas[1:]
                ind=fock_minus.index(new)
                coeff[ind]=psi0[i]
        
    if plus:
        coeff=utils.zeros(len(fock_plus))
        for i in xrange(len(fock)):
            bas=fock[i]
            if bas[0] == '0':
                new = '1'+bas[1:]
                ind=fock_plus.index(new)
                coeff[ind]=psi0[i]
    return coeff


def solve_perturb(nsites,U, hlat, omega,fock,e0,c0,sign_ham=1.):
    s=scipy.sparse.eye(len(fock))
    h=scipy.sparse.csr_matrix(MakeCorrelatedHam(len(fock),nsites,fock,U,hlat,hlat))
    lhs=omega*s-sign_ham*(h-e0*s)
    rhs=c0
    sol=scipy.sparse.linalg.spsolve(lhs,rhs)
    return sol

'''
def Evaluate_VPsi01(psi0,fock,plus,minus):
    fock_pert=[]        
    coeff_pert=[] 
    if minus:              
        for i in xrange(len(fock)):
            bas=fock[i]
            if bas[0] == '1':
                new = '0'+bas[1:]
                fock_pert.append(new)
                coeff_pert.append(psi0[i])        
    if plus:
        for i in xrange(len(fock)):
            bas=fock[i]
            if bas[0] == '0':
                new = '1'+bas[1:]
                fock_pert.append(new)
                coeff_pert.append(psi0[i])
    return N.array(coeff_pert), fock_pert
'''





def Swap(basis,j,k):
    new = basis[:j]+basis[k]+basis[j+1:k]+basis[j]+basis[k+1:]
    if k-j==1:
        phase=1.
    else:    
        phase = (-1.)**(sum_digits(int(basis[j+1:k])))
    return new, phase    
    
    


def MakeCorrelatedHam(l,L,fock,U,hup,hdw):
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
    return H_U+H_diag+H_t




def main():
    # main loop for single-particle GF
    utils.dtype=N.complex128
    nsites=8
    #t=1.
    cb_edge = 1.
    #rho0=(nsites-1)/(2*cb_edge)
    #gamma=0.1
    #t = N.sqrt(gamma/(rho0*N.pi))
    #gamma = (t**2)*N.pi*rho0
    #a=7.
    t = 0.1
    delta=0.005
    hlat=utils.zeros([nsites,nsites])
    hlat[0,1:] = t
    hlat[1:,0] = t
    cbstates = [0, 1./3, -1./3, 2./3, -2./3, 1., -1.]
    #cbstates = [-1., -2./3, -1./3, 0, 1./3, 2./3, 1.]
    for i in xrange(nsites-1):
        hlat[i+1,i+1] = cbstates[i]
    fock=GenerateFockSpace(nsites,False,False)
    l = len(fock)
    fock_minus=GenerateFockSpace(nsites,False,True)
    l_minus=len(fock_minus)
    fock_plus=GenerateFockSpace(nsites,True,False)
    l_plus=len(fock_plus)
    for U in [0.0, 2.5, 5.0, 7.5]:
	    fd=file("siam_dos_ed_nsites"+str(nsites)+"_U"+str(U)+"_eta"+str(delta),"w")
	    U*=t
	    hlat[0,0] = -0.5*U	  
	    Himp=scipy.sparse.csr_matrix(MakeCorrelatedHam(l,nsites,fock,U,hlat,hlat))
	    e,v=scipy.sparse.linalg.eigsh(Himp,k=1,which='SA')
	    E0=e[0]
	    psi0=v[:,0]
	    dpsi0=Evaluate_VPsi0(psi0,fock,fock_minus,fock_plus,False,True)
	    #dpsi0, fock_minus=Evaluate_VPsi01(psi0,fock,False,True)
	    cpsi0=Evaluate_VPsi0(psi0,fock,fock_minus,fock_plus,True,False)
	    #cpsi0, fock_plus=Evaluate_VPsi01(psi0,fock,True,False)
	    for omega in N.arange(-1.,1.,0.01):
		psi1=solve_perturb(nsites,U,hlat,omega+1j*delta,fock_minus,E0,dpsi0,sign_ham=-1.)
		G_d=N.dot(N.transpose(dpsi0.conjugate()),psi1)
		psi1=solve_perturb(nsites,U,hlat,omega+1j*delta,fock_plus,E0,cpsi0)
		G_c=N.dot(N.transpose(cpsi0.conjugate()),psi1)
		G=G_d+G_c
		print >>fd, omega, (-1./N.pi)*G.imag
        

#Calls the main function
if __name__ == "__main__":
    main()
       
