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
        # local time correlation function
        Ft = np.zeros(Ns)
 
        # Connectivity matirix for cluster
        #Cij = Mij.astype(int)  
        Cij = Mij.astype(int) # for multi-chain system
        Rank = np.zeros(1000) # Rank of the connectivity matrix
        clustsize = np.zeros(1000) # cluster size
        closedsticker = np.zeros(1000) # fraction of closed stickers
        opensticker = np.zeros(1000) # fraction of open stickers
        Pn = np.zeros(NAs) # binning to calculate probability dist of cluster size
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
        for tid, conf in enumerate(traj.box_coordinates):
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
                            Cij[tid, p, q] = 1
                q = 0
            if (tid % 1000 == 0):
                print 'tid = ', tid
          
        #computation of rank of connectivity matrix to find out number of cluster
        for t in range(Ns):
            if ((t+1) % 5 == 0):
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
                Rank[tid/5 - 1] = N
                clustsize[tid/5 - 1] = float(NAs)/N 
                for nc in range(N):
                    if (tid == 5):
                        print Cij[tid,:,nc].sum()
                    if (Cij[tid,:,nc].sum() > 1): 
                        closedsticker[tid/5 - 1] += Cij[tid,:,nc].sum()
                        if (Cij[tid,:,nc].sum() < NAs):
                            Pn[Cij[tid,:,nc].sum() - 1] += 1
                        else:
                            Pn[NAs-1] += 1
                    else:
                        opensticker[tid/5 - 1] += 1  
                        Pn[0] += 1
                #print Cij[tid,:,:]
                #print "#######################"
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
            np.savetxt("clustnum%s.txt" % filename, Rank)
            np.savetxt("clustsize%s.txt" % filename, clustsize)
            np.savetxt("closedstick%s.txt" % filename, closedsticker)
            np.savetxt("openstick%s.txt" % filename, opensticker)
            np.savetxt("clustdistb%s.txt" % filename, Pn)
        else:
            #print Ft
            #print gt 
            print Rank 
        
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
   for traj in range(1): 
       dir_name='/short/g16/as6030/traj_2.0/'
       #dir_name='/home/565/as6030/MixedFlow/MixedFlowCode1'
       base_filename=('r%sGamma0.1Chi1.0Nb24Nc21dt0.001Fraenkelb50.0phi2.0P5' % (traj+1))
       #base_filename='r0Gamma0.1Chi1.0Nb5Nc5dt0.001Fraenkelb50.0phi0.2P4'
       extension = 'nc'
       file = os.path.join(dir_name, base_filename + "." + extension)
       process(file, file_output = True, All_Prop = True, SF = False)

