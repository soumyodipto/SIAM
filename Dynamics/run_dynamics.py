import math
import time
import numpy as N
import scipy.linalg
import models
import models_embed
import embed
import qoperator
import response_embed
import opbasis
import opbasis_ni
import la
import utils

#import matplotlib.pyplot as plt

def main():
    # main loop for single-particle GF
    utils.dtype=N.complex128

    nimp=1
    nimp_sp=2*nimp	
    nsites=300
    nocc=nsites/2
    nsites_sp=nsites*2
    ndim=2
    sz=0
    # sites numbered in order as imp+bath ... ext
    sites_imp=range(0,2*nimp_sp)
    sites_ext=range(2*nimp_sp,nsites_sp)
    ncore=2*nocc-2*nimp
    cocc=range(2*nimp_sp,2*nimp_sp+ncore)
    vocc=range(2*nimp_sp+ncore,nsites_sp)

    N.set_printoptions(precision=3)

    #t=1. # hopping
    gamma=0.05
    rho0 = (nsites-2)/2.
    t = N.sqrt(gamma/(rho0*N.pi))
    #t = 0.1
    #gamma = (t**2)*pi*0.1
    delta=1e-12 # broadening
    #a=12.5	
#    for u in [0.0,4.0,10.0]:
    for u in [10.]:
        u *= gamma
        # Single impurity Anderson model
        mu=0.
        hlat,hall,hlat_sp,hall_sp,e0=models_embed.si_anderson_imp_ext_ni_h(t,u,nocc,nimp,nsites)
        hop,himp,hcs_imp,hds_imp,hcs_ext,hds_ext,hext=models_embed.si_anderson_imp_ext_h(hall_sp,u,nocc,nimp,nsites)
	#print hext
#        single site Hubbard in DMET basis
        #mu=u/2
        #hlat,hall,hlat_sp,hall_sp,e0=models_embed.hubbard_imp_ext_ni_h(t,u,nocc,nimp,nsites)
        #hop,himp,hcs_imp,hds_imp,hcs_ext,hds_ext,hext=models_embed.hubbard_imp_ext_h(hall_sp,u,nocc,nimp,nsites)

        # g.s embedding basis
        pemb,pcore,pvirt=embed.embed(hlat,nimp,nocc)

        # perturbation operator
        perturb=utils.zeros([nsites])
        perturb[0]=1.
        p_coeffs=N.array([1])
        perturb_dop=models.ContractedD(N.array([0]),p_coeffs)
        
        perturb_cop=models.ContractedC(N.array([0]),p_coeffs)

        #fd=file("siam_emb_dos_nsites"+str(nsites)+'_U'+str(u)+'_nimp'+str(nimp)+'_eta'+str(delta),"w")
        fdEnminus = file("En_minus","w")
	fdCnminus = file("Cn_minus","w")
	fdEnplus = file("En_plus","w")
	fdCnplus = file("Cn_plus","w")
        for omega in [1e-08]:#N.arange(-1.,1.,0.1):
    	    #print omega
            ops_dict=response_embed.mb_sp_ops(hlat,perturb,omega+1j*delta,nimp,nocc,pemb,pcore,pvirt)
            
            configs_dict=response_embed.mb_configs(nsites,nimp,nimp_sp,2*nocc-nimp_sp,0)
                      
            neutral_ops_configs,plus_ops_configs,minus_ops_configs=response_embed.mb_sp_ops_configs(ops_dict,configs_dict)
            
	    #print _ops_configs
            # basis is setup, now build matrix representations
            perturb_dop_mat=opbasis_ni.oimp_matrix_bra_ket_form(perturb_dop,minus_ops_configs,neutral_ops_configs,cocc,vocc)
            perturb_cop_mat=opbasis_ni.oimp_matrix_bra_ket_form(perturb_cop,plus_ops_configs,neutral_ops_configs,cocc,vocc)
            
            # h, N-1 configs
            
            #print ops_dict['d'][0]
            himp_mat_minus=opbasis_ni.oimp_matrix_form(himp,minus_ops_configs,cocc,vocc)            
            himp_ext_mat_minus=opbasis_ni.himp_ext_matrix_form(hcs_imp,hds_ext,hds_imp,hcs_ext,minus_ops_configs,cocc,vocc)
            #print himp_ext_mat_minus
            hext_mat_minus=opbasis_ni.hext_matrix_form(hext,minus_ops_configs,cocc,vocc,e0)
            hmat_minus=himp_mat_minus+himp_ext_mat_minus+hext_mat_minus
            unit_mat_minus=opbasis_ni.oimp_matrix_form(models.Unit(),minus_ops_configs,cocc,vocc)
            #print hext_mat_minus

            # h, N+1 configs
            himp_mat_plus=opbasis_ni.oimp_matrix_form(himp,plus_ops_configs,cocc,vocc)
            himp_ext_mat_plus=opbasis_ni.himp_ext_matrix_form(hcs_imp,hds_ext,hds_imp,hcs_ext,plus_ops_configs,cocc,vocc)
            hext_mat_plus=opbasis_ni.hext_matrix_form(hext,plus_ops_configs,cocc,vocc,e0)
            hmat_plus=himp_mat_plus+himp_ext_mat_plus+hext_mat_plus
            unit_mat_plus=opbasis_ni.oimp_matrix_form(models.Unit(),plus_ops_configs,cocc,vocc)

            # h, neutral configs
            himp_mat=opbasis_ni.oimp_matrix_form(himp,neutral_ops_configs,cocc,vocc)
            himp_ext_mat=opbasis_ni.himp_ext_matrix_form(hcs_imp,hds_ext,hds_imp,hcs_ext,neutral_ops_configs,cocc,vocc)
            hext_mat=opbasis_ni.hext_matrix_form(hext,neutral_ops_configs,cocc,vocc,e0)
            hmat=himp_mat+himp_ext_mat+hext_mat
            unit_mat=opbasis_ni.oimp_matrix_form(models.Unit(),neutral_ops_configs,cocc,vocc)
            
            #print hext_mat_minus
            
            # get neutral ground-state energy
            es,cs=la.peigh(hmat,unit_mat)
            e_gs=es[0]
            psi0=cs[:,0]
  	    
  	    
            es_minus,cs_minus=la.peigh(hmat_minus,unit_mat_minus)
            
    	    Enarray = es_minus - e_gs*N.ones(len(es_minus))
    	    print >>fdEnminus, omega, str(list(Enarray))[1:-1].translate(None, ',')
    	    es_plus,cs_plus=la.peigh(hmat_plus,unit_mat_plus)
    	    
    	    
    	    
    	    Enarray = es_plus - e_gs*N.ones(len(es_minus))
    	    print >>fdEnplus, omega, str(list(Enarray))[1:-1].translate(None, ',')

            gfminus=0.
    	    mat_el_list_minus=[]
    	    for i in xrange(len(es_minus)):
        	ket = N.dot(perturb_dop_mat,psi0)
        	bra = N.transpose(cs_minus[:,i].conjugate())
        	mat_el = N.dot(bra,ket)
        	mat_el *= mat_el.conjugate()
        	mat_el_list_minus.append(mat_el.real)
        	den = omega+1j*delta - e_gs + es_minus[i]
        	gfminus += (mat_el/den)

    	    print >>fdCnminus, omega, str(mat_el_list_minus)[1:-1].translate(None, ',')

    	    gfplus=0.
    	    mat_el_list_plus=[]
    	    for i in xrange(len(es_plus)):
		ket = N.dot(perturb_cop_mat,psi0)
		bra = N.transpose(cs_plus[:,i].conjugate())
		mat_el = N.dot(bra,ket)
		mat_el *= mat_el.conjugate()
		mat_el_list_plus.append(mat_el.real)
		den = omega+1j*delta + e_gs - es_plus[i]
		gfplus += (mat_el/den)

    	    print >>fdCnplus, omega, str(mat_el_list_plus)[1:-1].translate(None, ',')
            #print >>fd, omega, (-1./N.pi)*(gfplus.imag+gfminus.imag)
            #print (-1./N.pi)*(gfplus.imag+gfminus.imag)
            
                
#Calls the main function
if __name__ == "__main__":
    main()
