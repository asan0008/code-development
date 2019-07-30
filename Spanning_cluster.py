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
    print 'Lx = ', Lx
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

        pair = np.zeros((Ns, Nc*chains_ns[0])) # to define functionality in pairing

        spanning_matrix = np.zeros(500) # matrix to identify the system spanning network

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
                            # minimum image convention 
                            rijx = beadpositions[i][0] - beadpositions[j][0] - Lx*round((beadpositions[i][0] - beadpositions[j][0])/Lx)
                            rijy = beadpositions[i][1] - beadpositions[j][1] - Lx*round((beadpositions[i][1] - beadpositions[j][1])/Lx)
                            rijz = beadpositions[i][2] - beadpositions[j][2] - Lx*round((beadpositions[i][2] - beadpositions[j][2])/Lx)
                            rijsq = (rijx**2 + rijy**2 + rijz**2)
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
                 
                 
        # computation of the rank of the connectivity matrix for the chains and system spanning cluster
        for t, conf in enumerate(traj.box_coordinates):
            if ((t+1) % tdiff == 0):
                tid = t+1
                #print 'tid = ', tid
                for cid, chain_indices in enumerate(chains_indices):
                    positions = conf.array[chain_indices]
                    if cid == 0:
                        beadpositions = positions
                    else:
                        beadpositions = np.append(beadpositions, positions, axis=0) # writing position vectors of all the beads for all chains.
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
                    if (spanning_matrix[tid/tdiff - 1] == 0):                           
                        # identifying the chain number in a cluster
                        index = np.nonzero(Cij_chain[tid,:,nc])
                        index = index[0]
                        l = list(index) 
                        #print l
                        PBCx = 0; PBCy = 0; PBCz = 0
                        for j in range(len(l)):
                            for k in range(chains_ns[0]):
                                bd_i = chains_ns[0]*l[j] + k
                                if k == 0:
                                    xmin = beadpositions[bd_i][0]
                                    xmax = beadpositions[bd_i][0]
                                    ymin = beadpositions[bd_i][1]
                                    ymax = beadpositions[bd_i][1]
                                    zmin = beadpositions[bd_i][2]
                                    zmax = beadpositions[bd_i][2]
                                if k > 0:
                                    if (abs(beadpositions[bd_i][0] - beadpositions[bd_i-1][0]) > 50.0**0.5):
                                        if beadpositions[bd_i][0] > 0.:
                                            beadpositions[bd_i][0] = -Lx + beadpositions[bd_i][0]
                                        elif beadpositions[bd_i][0] < 0.:
                                            beadpositions[bd_i][0] = Lx + beadpositions[bd_i][0]
                                    
                                    if (abs(beadpositions[bd_i][1] - beadpositions[bd_i-1][1]) > 50.0**0.5):
                                        if beadpositions[bd_i][1] > 0.:
                                            beadpositions[bd_i][1] = -Lx + beadpositions[bd_i][1]
                                        elif beadpositions[bd_i][1] < 0.:
                                            beadpositions[bd_i][1] = Lx + beadpositions[bd_i][1]
                                    
                                    if (abs(beadpositions[bd_i][2] - beadpositions[bd_i-1][2]) > 50.0**0.5):
                                        if beadpositions[bd_i][2] > 0.:
                                            beadpositions[bd_i][2] = -Lx + beadpositions[bd_i][2]
                                        elif beadpositions[bd_i][2] < 0.:
                                            beadpositions[bd_i][2] = Lx + beadpositions[bd_i][2]
                                    # calculation max and min position co-ordinates of a chain
                                    if beadpositions[bd_i][0] < xmin: 
                                        xmin = beadpositions[bd_i][0]
                                    elif beadpositions[bd_i][0] > xmax:
                                        xmax = beadpositions[bd_i][0]
                                      
                                    if beadpositions[bd_i][1] < ymin:
                                        ymin = beadpositions[bd_i][1]
                                    elif beadpositions[bd_i][1] > ymax:
                                        ymax = beadpositions[bd_i][1]

                                    if beadpositions[bd_i][2] < zmin:
                                        zmin = beadpositions[bd_i][2]
                                    elif beadpositions[bd_i][2] > zmax:
                                        zmax = beadpositions[bd_i][2]
                                #print 'k = ', k
                                #print 'xmin = ', xmin
                                #print 'xmax = ', xmax
                             
                            if xmin < -0.5*Lx:
                                xmin1 = xmin; xmin2 = xmin + Lx; xmax1 = xmax; xmax2 = xmax + Lx
                                PBCx = 1
                            if xmax > 0.5*Lx:
                                xmin1 = xmin - Lx; xmin2 = xmin; xmax1 = xmax - Lx; xmax2 = xmax
                                PBCx = 1
                            elif xmin > -0.5*Lx and xmax < 0.5*Lx:
                                xmin1 = xmin - Lx; xmin2 = xmin; xmax1 = xmax; xmax2 = xmax + Lx                               
 
                            if ymin < -0.5*Lx:
                                ymin1 = ymin; ymin2 = ymin + Lx; ymax1 = ymax; ymax2 = ymax + Lx
                                PBCy = 1
                            if ymax > 0.5*Lx:
                                ymin1 = ymin - Lx; ymin2 = ymin; ymax1 = ymax - Lx; ymax2 = ymax
                                PBCy = 1
                            elif ymin > -0.5*Lx and ymax < 0.5*Lx:
                                ymin1 = ymin - Lx; ymin2 = ymin; ymax1 = ymax; ymax2 = ymax + Lx

                            if zmin < -0.5*Lx:
                                zmin1 = zmin; zmin2 = zmin + Lx; zmax1 = zmax; zmax2 = zmax + Lx
                                PBCz = 1
                            if zmax > 0.5*Lx:
                                zmin1 = zmin - Lx; zmin2 = zmin; zmax1 = zmax - Lx; zmax2 = ymax
                                PBCz = 1
                            elif zmin > -0.5*Lx and zmax < 0.5*Lx:
                                zmin1 = zmin - Lx; zmin2 = zmin; zmax1 = zmax; zmax2 = zmax + Lx  
                            
                            if (j == 0):                          
                                xmin1c = xmin1; xmin2c = xmin2; xmax1c = xmax1; xmax2c = xmax2
                                ymin1c = ymin1; ymin2c = ymin2; ymax1c = ymax1; ymax2c = ymax2
                                zmin1c = zmin1; zmin2c = zmin2; zmax1c = zmax1; zmax2c = zmax2
                            
                            if j > 0:
                                xmax2c = min(xmax2c, xmax2); xmin1c = max(xmin1c, xmin1)
                                xmax1c = max(xmax1c, xmax1); xmin2c = min(xmin2c, xmin2)

                                ymax2c = min(ymax2c, ymax2); ymin1c = max(ymin1c, ymin1)
                                ymax1c = max(ymax1c, ymax1); ymin2c = min(ymin2c, ymin2)                             
                               
                                zmax2c = min(zmax2c, zmax2); zmin1c = max(zmin1c, zmin1)
                                zmax1c = max(zmax1c, zmax1); zmin2c = min(zmin2c, zmin2)
                            '''
                            print 'j = ', j
                            print 'xmin1c = ', xmin1c
                            print 'xmin2c = ', xmin2c
                            print 'xmax1c = ', xmax1c
                            print 'xmax2c = ', xmax2c
                            '''
                            if xmax1c >= xmin2c or ymax1c >= ymin2c or zmax1c >= zmin2c:
                                if PBCx == 0:
                                    clustrx = abs(xmax1c - xmin2c)
                                elif PBCx == 1:
                                    clustrx = min(abs(xmax2c - xmin2c), abs(xmax1c - xmin1c))
                                if PBCy == 0:
                                    clustry = abs(ymax1c - ymin2c)
                                elif PBCy == 1:
                                    clustry = min(abs(ymax2c - ymin2c), abs(ymax1c - ymin1c))
                                if PBCz == 0:
                                    clustrz = abs(zmax1c - zmin2c)
                                elif PBCz == 1:
                                    clustrz = min(abs(zmax2c - zmin2c), abs(zmax1c - zmin1c))
 
                                if clustrx >= Lx-1 or clustry >= Lx-1 or clustrz >= Lx-1:
                                    #print 'clustrx = ', clustrx
                                    #print 'clustry = ', clustry
                                    #print 'clustrz = ', clustrz
                                    spanning_matrix[tid/tdiff - 1] = 1
                                    break
                             
        if file_output: 
            np.savetxt("spanningclust%s.txt" % filename, spanning_matrix)
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
   for traj in range(25,48): 
       dir_name='/short/g16/as6030/traj_f4cbycstr0.3phi5.0/'
       #dir_name='/home/565/as6030/MixedFlow/MixedFlowCode1'
       base_filename=('r%sGamma0.1Chi1.0Nb24Nc33dt0.001FENEb50.0phi5.0conc0.3' % (traj+1))
       #base_filename='r0Gamma0.1Chi1.0Nb5Nc5dt0.001Fraenkelb50.0phi0.2P4'
       extension = 'nc'
       file = os.path.join(dir_name, base_filename + "." + extension)
       process(file, file_output = True, All_Prop = True, SF = False)

