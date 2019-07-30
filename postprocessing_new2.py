#!/usr/bin/env python
import os

import numpy as np

from MMTK import *

from MMTK.Trajectory import Trajectory
from MMTK import AtomCluster

execfile ('Input.py')
def process(cdffile, file_output, All_Prop, SF):
    dirname, filepath = os.path.split(cdffile)
    filename, ext = os.path.splitext(filepath)

    traj = Trajectory(None, cdffile)
    universe = traj.universe
    forces = traj.gradients
    box_posvec = traj.box_coordinates
    chains = universe.objectList(AtomCluster)
    chains_indices = [[atom.index for atom in chain.atomList()] for chain in chains]
    chains_ns = [chain.numberOfAtoms() for chain in chains]
    print ('In postprocessing %s' % filename)
    print 'ASPeriodicity = ', periodicity
    print 'Lja2 = ', Lja2
    # Number of sample points
    Ns = len(traj)
    # Number of chains
    Nc = len(chains)
    print (Ns, Nc, chains_ns[0])
    #NAs = (int(chains_ns[0]/periodicity) + 1)*Nc # total no. of associating beads for all the chains for telechelic system
    NAs = (int(chains_ns[0]/periodicity))*Nc # total no. of associating beads for all the chains for multi-sticker system
    print NAs
    tdiff = 4 # N_prdt/tdiff : no. of data points for the postprocessing 

    if All_Prop:
        # center of mass
        r_coms = np.zeros((Ns, Nc, 3))
        # end to end distance
        q_sqs = np.zeros((Ns, Nc))
        # Radius of gyration
        rg_sqs = np.zeros((Ns, Nc))
        # tau xy
        tauxys = np.zeros(Ns)

        ## Neighbour matrix Mij(for association dynamics), added by Aritra
        #Mij = np.zeros((Ns, chains_ns[0], chains_ns[0])) # for single chain system
        Mij = np.zeros((Ns, NAs, NAs)) # for multi-chain system
        Mij_chain = np.zeros((Ns, Nc, Nc)) # to calculate connectivity matrix for chains
        # local time correlation function
        Ft = np.zeros(Ns)
 
        # Connectivity matirix for cluster
        #Cij = Mij.astype(int)  
        Cij = Mij.astype(int) # for multi-chain system
        Cij_chain = Mij_chain.astype(int) 

        Rank = np.zeros(500) # Rank of the connectivity matrix
        Rank_chain = np.zeros(500) # Rank of the connectivity matrix for chains
        clustsize = np.zeros(500) # cluster size

        closedsticker = np.zeros(500) # number of closed stickers
        bondedchain_inter = np.zeros(500) # number of inter-chain cluster
        opensticker = np.zeros(500) # number of open stickers
        bondedchain_intra = np.zeros(500) # number of intra-chain cluster

        Pn = np.zeros(NAs) # binning to calculate probability dist of cluster size
        Pn_chain = np.zeros(Nc) # binning to calculate probability dist of chain cluster size

        open_beads = np.zeros(Nc) # to calculate open beads per chain
        intrachain_beads = np.zeros(Nc) # to calculate intrachain associated beads per chain
        openchain = np.zeros(500) # to calculate number of  open chains
        interchain = np.zeros(500) # to calculate number of interchain association
        pair = np.zeros((Ns, Nc*chains_ns[0])) # to define functionality in pairing
        ''' 
        for tid, conf in enumerate(traj.configuration):
            tauxy = 0.0
            for cid, chain_indices in enumerate(chains_indices):
                positions = conf.array[chain_indices]
                r_com = positions.sum(axis = 0) / chains_ns[cid]
                q_sq = ((positions[0] - positions[-1])**2).sum()
                rel_pos = positions - r_com
                rg_sq = (rel_pos**2).sum() / chains_ns[cid]
                tauxy += (rel_pos[:, 0] * forces[tid].array[chain_indices][:, 1]).sum()
                r_coms[tid, cid] = r_com
                q_sqs[tid, cid] = q_sq
                rg_sqs[tid, cid] = rg_sq
            tauxys[tid] = tauxy

        # computation of Neighbour matrix for single chain
        for tid, conf in enumerate(traj.configuration):
            for cid, chain_indices in enumerate(chains_indices):
                positions = conf.array[chain_indices]
                for i in range(chains_ns[cid]-2):
                    for j in range(i+2, chains_ns[cid]):
                        rijsq = ((positions[i] - positions[j])**2).sum()
                        rij = sqrt(rijsq)
                        if rij < 1.5:
                           Mij[tid, i, j] = 1
        '''
        ''' 
        # computation of Neighbour matrix for multi-chain systems
        for tid, conf in enumerate(traj.box_coordinates):
            for cid, chain_indices in enumerate(chains_indices):
                positions = conf.array[chain_indices]
                if cid == 0:
                    beadpositions = positions
                else:
                    beadpositions = np.append(beadpositions, positions, axis=0) # writing position vectors of all the beads for all chains. 
            p = 0 # bead index for associating bead i
            q = 0 # bead index for associating bead j
            for i in range(chains_ns[0]*Nc - 1):
                nc_i = int(i/chains_ns[0]) # chain number for bead i
                if ((i - nc_i) % periodicity == 0 and i > 0):
                    p += 1
                    if (p > (NAs-1)):
                        print 'p is greater than NAs: ', p+1
                for j in range(i+1, chains_ns[0]*Nc):
                    nc_j = int(j/chains_ns[0]) # chain number for bead j 
                    if ((j - nc_j) % periodicity == 0 and j > 0):
                        q = p + 1
                        if (q > (NAs-1)):
                            print 'q is greater than NAs: ', q+1
                    if (((i - nc_i) % periodicity) == 0 and ((j - nc_j) % periodicity) == 0):
                        rijsq = ((beadpositions[i] - beadpositions[j])**2).sum()
                        rij = sqrt(rijsq)
                        if (rij < 1.5*Lja2):
                            Mij[tid, p, q] = 1
                q = 0
            if (tid % 2000 == 0):
                print 'tid = ', tid    
  
        ## compution of number of cluster and cluster size for single chain system
        for tid, conf in enumerate(traj.configuration):
            for cid, chain_indices in enumerate(chains_indices):
                positions = conf.array[chain_indices]
                for i in range(chains_ns[cid]):
                    for j in range(chains_ns[cid]):
                        rijsq = ((positions[i] - positions[j])**2).sum()
                        rij = sqrt(rijsq)
                        if rij < 2.8:
                           Cij[tid, i, j] = 1  
        # computation of rank of connectivity matrix to find out number of cluster
        for t in range(Ns):
            if ((t+1) % 1000 == 0):
                tid = t+1 
                N0 = chains_ns[0]
                N = chains_ns[0]
                #print Cij[tid,:,:]
                #print "****************"
                for i in range(N0):
                    if i < N:
                        switch = 1
                        while (switch == 1):
                            for j in range(i, N):
                                switch = 0
                                c = np.bitwise_and(Cij[tid,:,i], Cij[tid,:,j]) 
                                if ((c**2).sum() != 0 and j!=i):
                                    switch = 1
                                    Cij[tid,:,i] =  np.bitwise_or(Cij[tid,:,i], Cij[tid,:,j]) 
                                    N1 = j
                                    N = N-1
                                    for k in range(N1+1, N+1):
                                        Cij[tid,:,k-1] = Cij[tid,:,k]
                                    break     
                    else: break 
                Rank[tid/1000 - 1] = N
                clustsize[tid/1000 - 1] = float(chains_ns[0])/N
                #print Cij[tid,:,:]
                #print "#######################"
        '''
        ## computation of number of cluster and cluster size for multi-chain system 
        for t, conf in enumerate(traj.box_coordinates):
            if ((t+1) % tdiff == 0):
                tid = t+1
                #print tid
                for cid, chain_indices in enumerate(chains_indices):
                    positions = conf.array[chain_indices]
                    if cid == 0:
                        beadpositions = positions
                    else:
                        beadpositions = np.append(beadpositions, positions, axis=0) # writing position vectors of all the beads for all chains.
                p = 0 # bead index for associating bead i
                q = 0 # bead index for associating bead j
                for i in range(chains_ns[0]*Nc):
                    nc_i = int(i/chains_ns[0]) # chain number for bead i
                    # for telechelic system
                    '''
                    if ((i - nc_i) % periodicity == 0 and i > 0):
                        p += 1
                        if (p > (NAs-1)):
                            print 'p is greater than NAs: ', p+1
                    '''
                    # for multi-sticker system (4 : periodicity-1)
                    if (((i-4)+nc_i) % periodicity == 0 and i > 4):
                        p += 1
                        if (p > (NAs-1)):
                            print 'p is greater than NAs: ', p+1
                    for j in range(chains_ns[0]*Nc):
                        nc_j = int(j/chains_ns[0]) # chain number for bead j
                        # for telechelic system
                        '''
                        if ((j - nc_j) % periodicity == 0 and j > 0):
                            q += 1
                            if (q > (NAs-1)):
                                print 'q is greater than NAs: ', q+1
                        '''
                        # for multi-sticker system (4: periodicity-1)
                        if (((j-4)+nc_j) % periodicity == 0 and j > 4):
                            q += 1
                            if (q > (NAs-1)):
                                print 'q is greater than NAs: ', q+1
                        # for telechelic system
                        '''
                        if (((i - nc_i) % periodicity) == 0 and ((j - nc_j) % periodicity) == 0):
                            rijsq = ((beadpositions[i] - beadpositions[j])**2).sum()
                            rij = sqrt(rijsq)
                            if (rij < 1.5*Lja2):
                                Cij[tid, p, q] = 1
                        '''
                        # for multi-sticker system
                        if ((((i-4)+nc_i) % periodicity) == 0 and (((j-4)+nc_j) % periodicity) == 0):
                            rijsq = ((beadpositions[i] - beadpositions[j])**2).sum()
                            rij = sqrt(rijsq)
                            if (rij < 1.82*Lja2):
                                if (i == j):
                                    Cij[tid, p, q] = 1
                                    Cij_chain[tid, nc_i, nc_j] = 1
                                if (i != j and pair[tid, i] < 1 and pair[tid, j] < 1):
                                    Cij[tid, p, q] = 1
                                    Cij_chain[tid, nc_i, nc_j] = 1
                                    pair[tid, i] += 1
                                    pair[tid, j] += 1
                    q = 0
                if (tid % 500 == 0):
                    print 'tid = ', tid
        '''          
        #computation of rank of connectivity matrix to find out number of cluster
        for t in range(Ns):
            if ((t+1) % tdiff == 0):
                tid = t+1
                N0 = NAs
                N = NAs
                #print Cij[tid,:,:]
                #print "****************"
                for i in range(N0):
                    if i < N:
                        switch = 1
                        while (switch == 1):
                            for j in range(i, N):
                                switch = 0
                                c = np.bitwise_and(Cij[tid,:,i], Cij[tid,:,j])
                                if ((c**2).sum() != 0 and j!=i):
                                    switch = 1
                                    Cij[tid,:,i] =  np.bitwise_or(Cij[tid,:,i], Cij[tid,:,j])
                                    N1 = j
                                    N = N-1
                                    for k in range(N1+1, N+1):
                                        Cij[tid,:,k-1] = Cij[tid,:,k]
                                    break
                    else: break
                Rank[tid/tdiff - 1] = N
                clustsize[tid/tdiff - 1] = float(NAs)/N 
                for nc in range(N):
                    if (Cij[tid,:,nc].sum() > 1):
                        
                        #if tid == 20:
                         #   print "col no. = ", nc
                         #   print Cij[tid,:,nc]
                            
                        closedsticker[tid/tdiff - 1] += Cij[tid,:,nc].sum()
                        if (Cij[tid,:,nc].sum() < NAs):
                            Pn[Cij[tid,:,nc].sum() - 1] += 1
                        else:
                            Pn[NAs-1] += 1
                        flag1 = 0
                        flag2 = 0 
                        # this part calculates the fraction of beads which are inter-chain associated
                        for n in range(Cij[tid,:,nc].size):
                            if (Cij[tid,n,nc] != 0 and flag2 == 0):
                                chain_num0 = int(n/(NAs/Nc))
                                flag1 += 1
                                if flag1 == 1:
                                   chain_num1 = chain_num0
                                if (flag1 > 1 and chain_num0 != chain_num1):
                                   interchain[tid/tdiff - 1] += Cij[tid,:,nc].sum()
                                   flag2 = 1
                            elif (flag2 == 1): break 
                        if (flag2 == 0):
                            intrachain_beads[chain_num1] += Cij[tid,:,nc].sum()                       
                            
                    else:
                        opensticker[tid/tdiff - 1] += 1  
                        Pn[0] += 1
                        # calculating fraction of open chains for multi-sticker polymers
                        for n in range(Cij[tid,:,nc].size):
                            if Cij[tid,n,nc] != 0:    
                                chain_num = int(n/(NAs/Nc)) 
                                open_beads[chain_num] += 1
                #print Cij[tid,:,:]
                #print "#######################"
                for cn in range(Nc):
                    #print tid
                    #print open_beads[cn]
                    intrachain_beads[cn] += open_beads[cn]
                    
                    # if (tid == tdiff):
                        # print "intrachain_%s = %s" % (cn, intrachain_beads[cn]) 
                    
                    if (open_beads[cn] == NAs/Nc):
                        openchain[tid/tdiff - 1] += 1 
                    open_beads[cn] = 0
                    intrachain_beads[cn] = 0
        '''         
        # computation of the rank of the connectivity matrix for the chains
        for t in range(Ns):
            if ((t+1) % tdiff == 0):
                tid = t+1
                N0 = Nc
                N = Nc
                #if (tid == 8):
                #    print Cij_chain[tid,:,:]
                #    print "****************"
                for i in range(N0):
                    if i < N:
                        switch = 1
                        while (switch == 1):
                            for j in range(i, N):
                                switch = 0
                                c = np.bitwise_and(Cij_chain[tid,:,i], Cij_chain[tid,:,j])
                                if ((c**2).sum() != 0 and j!=i):
                                    switch = 1
                                    Cij_chain[tid,:,i] =  np.bitwise_or(Cij_chain[tid,:,i], Cij_chain[tid,:,j])
                                    N1 = j
                                    N = N-1
                                    for k in range(N1+1, N+1):
                                        Cij_chain[tid,:,k-1] = Cij_chain[tid,:,k]
                                    break
                    else: break
                Rank_chain[tid/tdiff - 1] = N
                #if (tid == 8):
                #    print Cij_chain[tid,:,:]
                for nc in range(N):
                    if (Cij_chain[tid,:,nc].sum() > 1):
                        '''
                        if tid == 20:
                            print "col no. = ", nc
                            print Cij[tid,:,nc]
                        '''
                        bondedchain_inter[tid/tdiff - 1] += Cij_chain[tid,:,nc].sum()
                        if (Cij_chain[tid,:,nc].sum() < Nc):
                            Pn_chain[Cij_chain[tid,:,nc].sum() - 1] += 1
                        else:
                            Pn_chain[Nc-1] += 1
                        
                    else:
                        bondedchain_intra[tid/tdiff - 1] += 1
                        Pn_chain[0] += 1
                 
        '''                       
        # mean square displacement
        msds = np.zeros((Ns, Nc))
        for dtime in range(1, Ns):
            for cid in range(Nc): 
                com = r_coms[:, cid]
                msd = ((com[dtime:, :] - com[:-dtime, :])**2).sum()
                msd /= Ns - dtime 
                msds[dtime, cid] = msd
       
        # stress autocorrelation 
        gt = np.zeros(Ns)
        gt[0] = (tauxys**2).sum() / (Nc * Ns)
        for dtime in range(1, Ns):
            gt[dtime] = (tauxys[dtime:] * tauxys[:-dtime]).sum() / (Nc * (Ns - dtime))
            if (dtime % 2000 == 0):
                print 'dtime = ', dtime
        
        # local time correlation function for association dynamics
        Ft[0] = (Mij**2).sum()/(Nc*Ns*chains_ns[0]) 
        print 'dtime = 0' 
        for dtime in range(1, Ns):
            Ft[dtime] = (Mij[dtime:]*Mij[:-dtime]).sum()/(Nc*(Ns - dtime)*chains_ns[0])
            if (dtime % 2000 == 0):
                print 'dtime = ', dtime
        '''
        #data_out = np.column_stack( (traj.time, np.mean(q_sqs, axis = 1), np.mean(rg_sqs, axis = 1), np.mean(msds, axis = 1), gt, tauxys, Ft) )
        
        '''
        Rg2 = np.mean(rg_sqs, axis = 1)
        Rg2mean = np.mean(Rg2)
        #print Rg2mean
        with open('Rg2mean.txt','a+') as fRg2:
              fRg2.write("%lf\n" % Rg2mean)
        '''
        if file_output:
            #np.savetxt("dynamic%s.txt" % filename, data_out)
            #np.savetxt("Rg2.txt", Rg2)
            #np.savetxt("Ft%s.txt" % filename, Ft)
            #np.savetxt("Gt%s.txt" % filename, gt)
            #np.savetxt("clustnum%s.txt" % filename, Rank)
            np.savetxt("chainclustnum%s.txt" % filename, Rank_chain)
            #np.savetxt("closedstick%s.txt" % filename, closedsticker)
            np.savetxt("bondedchain%s.txt" % filename, bondedchain_inter)
            #np.savetxt("openstick%s.txt" % filename, opensticker)
            np.savetxt("unbondedchain%s.txt" % filename, bondedchain_intra)
            #np.savetxt("clustdistb%s.txt" % filename, Pn)
            np.savetxt("chainclustdistb%s.txt" % filename, Pn_chain)   
            #np.savetxt("openchain%s.txt" % filename, openchain)
            #np.savetxt("freechain%s.txt" % filename, openchain)
            #np.savetxt("interchain%s.txt" % filename, interchain)  
        else:
            #print Ft
            #print gt 
            print Rank_chain 
        
    if SF:
        # structure factor
        from math import sin
        structure_kmin, structure_kmax, structure_nks = 0.1, 8, 100
        ks = np.logspace(np.log10(structure_kmin), np.log10(structure_kmax), structure_nks)
        structure_factors = np.zeros((structure_nks, Ns, Nc))
        Npair = ((Nbpc*Nbpc) - Nbpc)/2
        r_mag = np.zeros((Ns, Nc, Npair))
        Ns = len(traj)
        for tid, conf in enumerate(traj.configuration):
            for cid, chain_indices in enumerate(chains_indices):
                positions = conf.array[chain_indices]
                n = len(positions)
		pairid = 0
                for i, position1 in enumerate(positions):
                    for j in range(i+1, n):
                        position2 = positions[j]
                        rij = position2 - position1
                        r_mag[tid, cid, pairid] = ((rij * rij).sum())**.5
			pairid = pairid + 1
#########################################################################
        for kid, k in enumerate(ks):
            print 'kid', kid
            for tid in range(Ns):
                for cid in range(Nc):
                    struct_sum = 0.0
		    pairid = 0
                    for i in range(Nbpc):
                        for j in range(i+1, Nbpc):
                            struct_sum += sin(k*r_mag[tid, cid, pairid])/(k*r_mag[tid, cid, pairid])
		            pairid = pairid + 1
                    structure_factors[kid, tid, cid] = 1.0 + (2 * struct_sum / Nbpc) #multiply by two because of symmetricity and addition of 1 is because of rij=0 terms
        print 'kid', kid
        structure_factors = np.mean(structure_factors, axis = 2)
        data_structure_factor = np.column_stack( (ks, np.mean(structure_factors, axis = 1), np.std(structure_factors, ddof = 1, axis = 1)) )
        if file_output:
            np.savetxt("structure%s.txt" % filename, data_structure_factor)

    traj.close()

if __name__ == "__main__":
   # import glob
   # files = glob.glob("r*.nc")
   # for file in files:
   for traj in range(50): 
       dir_name='/short/g16/as6030/traj_f5cbycstr0.3/'
       #dir_name='/home/565/as6030/MixedFlow/MixedFlowCode1'
       base_filename=('r%sGamma0.1Chi1.0Nb29Nc9dt0.001FENEb50.0phi5.0P5' % (traj+1))
       #base_filename='r0Gamma0.1Chi1.0Nb5Nc5dt0.001Fraenkelb50.0phi0.2P4'
       extension = 'nc'
       file = os.path.join(dir_name, base_filename + "." + extension)
       process(file, file_output = True, All_Prop = True, SF = False)

