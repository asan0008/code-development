#from math import *
from math import sin, cos, pi, sqrt, log, atan, log10
import random
#################################################################################################
## There are 6 sections in this input file. In each section user needs to specify the following:#
## Section 1. Supply the values of some flags according to usage				#
## Section 2. Supply input parameters related to the system, polymer solution etc..		#
## Section 3. Supply flow parameters								#
## Section 4. Supply parameters related to all the interactions/forces				#
## Section 5. Supply parameters related to job length, time steps, statistics etc..		#
## Section 6: Specify file names for the input/output						#
#################################################################################################
## SECTION 1: SOME FLAGS									#
#################################################################################################
#BATCH = True # FOR JOB SUBMISSIONS
BATCH = False # FOR DIRECT RUNS/for equilibriation
#SERIES = True
#SERIES = False #default to False (make True if you want run long jobs "in series" - e.g.-NCI runs)
#APPEND = False #default to False (do not change unless you want to continue a jobs - requires you to specify the file to be appended) 
#APPEND = True # if appending is done in series
#SINGLE=True 
#SINGLE=False
#PREEQUILIBRATION = True # if preequilibration is true: make new equilibrated chain from scratch
PREEQUILIBRATION = False  # otherwise: use already equilibrated chain as starting configuration
file_output = True # creates a .txt output file
#file_output = False
All_Prop = True #calculates all propertise (excludes structure factor which has to be specified below)
#All_Prop = False
#SF = True  #calculates the structure factor as well
SF = False
#######################
## FLOW FLAGS
## Only Shear Flow: SHEAR = 1 and ELONGATION = 0
## Only Elongational Flow: SHEAR = 0 and ELONGATION = 1
## Mixed Flow: SHEAR = 1 and ELONGATION = 1
SHEAR = 0
ELONGATION = 0 
EQinEQ = 1 # 1 if equilibration is to be done without flow and 0 if flow has to happen
#######################
# Rejection Algorithm (RA) Flag
RA = 1 # 0 if no RA and 1 if RA
# Rejection Algorithm Parameters 
timeStepFactor = 1.0 # A factor by which the time step is reduced or ramped up
N_wait_rampup = 50 # Number of times to be waiting for the time step to be ramped up
##################################################################################################
## SECTION 2: POLYMER SOLUTION/SYTEM PARAMETERS							 #
##################################################################################################
cbycstar = 0.3  # c/c*
Nbpc = 29 # Number of beads per chain
a = 0.322738 # Non-dimensional bead radius (a = dim. bead rad/l_H)
springType = 1 #0: Hookean, 1:FENE, 2: ILC, 3:WLC, 4:FENE30 (Steinhauser potential), 5:Fraenkel
b_para = 50.0  # FENE b parameter OR |Q0|^2
Q_para = 1.0 # natural length for Fraenkel spring
#####################
## box size and number of chains
#Rg0 = sqrt(((Nbpc*Nbpc) - 1.0)/(2.0*Nbpc)) #valid only for theta solvent (no EV) radius of gyration
#Re0 = sqrt(3.0*(Nbpc - 1.0)) #valid only for theta solvent (no EV) end to end distance
Rg0 = sqrt(12.726) 
Re0 = sqrt(76.356)
Rg = Rg0
exponent = ((2.0*0.6) - 1.0) / (1.0 - (3.0*0.6))
Re = sqrt(Re0*Re0*(cbycstar**exponent))
#L = 2.0*Re
L = 2.0*Re0
Nc = int(cbycstar*3.0*L*L*L/(4.0*pi*Rg*Rg*Rg))+1
#Nc = 10
# Initial Box Sizes (Starts with Orthogonal system)
Lx = L 
Ly = 1.0*L
Lz = 1.0*L
#################################################################################################
# SECTION 3: Flow Parameters									#
#################################################################################################	
Chi = 1.0  # Mixedness parameter
Gamma = 0.1 
k_para = 3 # No need to change
lamdap_MF = (k_para + (sqrt((k_para*k_para) - 4)))/2 # No need to change
eps_p_MF = log(lamdap_MF) # No need to change
N11 = 2 # No need to change
N12 = -sqrt((N11*(k_para-N11)) - 1) # No need to change
theta_MF = atan((N11 - lamdap_MF)/N12) # Inclination of the box
gam_dot = (1-Chi)*Gamma # Shear rate
if ELONGATION == 1:
	eps_dot = Gamma*sqrt(Chi) # Elongational rate
	taup_MF = eps_p_MF/eps_dot # Strain period for extensional flow
else:
	eps_dot = 0.0000000001 # Dummy and irrelevant for shear flow
        taup_MF = eps_p_MF/eps_dot # Dummy and irrelevant for shear flow
taup_SH = 0.0 #Dummy value for taup_SH (Irrelevant for extensional or mixed flow)
###############################################################################################
## SECTION 4: Parameters for all the interactions (HI, EV, DH, ES) along with their flags     #
###############################################################################################
# Non-hydrodynamic Forces 
#############################
## LJ/Narrow Gaussian PARAMETERS
######################
typeEV = 0 ## 2 for narrow Gaussian, 1 for capped Lenard Jones, 3 for full Lennard Jones, and 0 for NO EV
subbox_per_rcEV = 2 # Neighborlist parameter
## LJ
kfene_bar = 7.0 # Non-dimensional kfene in LB units (Tri et al JCP 2009)
ktbyepsilon = 10.0 # Non-dimensional temperature in LB units (Tri et al JCP 2009)
sigmabylk = sqrt(kfene_bar/ktbyepsilon) # Eq 8 (Tri et al JCP 2009)
#ktbyepsilon = 1.2 # for stoltz
#sigmabylk = 1
lj1 = ktbyepsilon
lj2 = sigmabylk
#Dstar = lj2
#Zstar = 4.0/lj1
#Hstar = a/sqrt(pi)
#EV_cutoff = L / 5. # Eq 2 (Tri et al JCP 2009)
#######################
## Narrow Gaussian
#######################
z = 1.0
#zst=z/sqrt(Nbpc)  ##enter the value of dimensionless variable z* here
zst = 0.440417
dst=zst**0.2
#dst = 1
EV_cutoff = 5*dst
#EV_cutoff = L/10.0

#############################
## ASSOCIATION   (Incorporated on Aug 2017 by Aritra)
######################
typeAS = 1  #1 for Modified Lennard Jones, 0 for NO ASSOCIATION
periodicity = 5 # indicates the position of the associating beads along the backbone of the chain
subbox_per_rcAS = 2 # Neighborlist parameter
## LJ interation between associative groups
Lja1 = 1.0 # LJ parameter for association
Lja2 = 1.0 # LJ parameter for associaiton
phi = 5.0 # depth of the well
alphaAS = 1.530633312 # parameter determined by T. Soddemann, B. Duenweg, 2001 with cutoff 1.82
betaAS = 1.213115524 # parameter determined by T. Soddemann, B. Duenweg, 2001 with cutoff 1.82 
'''
alphaAS = 3.1730728678 # parameter determined by T. Soddemann, B. Duenweg, 2001 with cutoff 1.5
betaAS = -0.85622864544 # parameter determined by T. Soddemann, B. Duenweg, 2001 with cutoff 1.5
'''
AS_cutoff1 = (2.**(1./6.))*Lja2 # end of repulsive domain  
AS_cutoff2 = 1.82*Lja2  # end of attractive domain

#######################
## charges
######################
chargename = 'bd_charge'    # name of the charge property
charge = 0.0               # charge on beads in e
#######################
## Debye-Huckel (DH) PARAMETER (dimensionless)
######################
typeDH = 0           # 1 for having debye huckel and 0 for no debye
subbox_per_rcDH = 1  # Neighborlist parameter. This number is not recommended to change.
Lbstar=0.46          # Lb*=(lb/bk) * bk*  where bk*=bk/lh where lh is the scaling factor
Ldstar=26.12         # Ld*=(ld/bk) * bk*
#dncharge=      #charge on a DNA bead
#prcharge=      #charge on a protein bead
DH_cutoff=0.98*L      # cutoff radius for Debye
lowerlimits=lj2*DH_cutoff
#######################
## Electrostatics (ES) parameter
######################
typeES = 0              # 0 for no ES, 1 for ES
subbox_per_rcES = 1     # Neighborlist parameter. This number is not recommended to change.
energy_factorES = 0  # factor for the electrostatic energy responsible for unit conversions
######################
## Hydrodynamic interactions
######################
typeHI = 2 # 1 for Hydrodynamic interaction 2 for without hydrodynamic interaction
## Ewald parameter
M = 3.0 # Accuracy parameter
factor1 = 0.25 #Ewald optimization parameter
factor2 = 0.1 #Ewald optimization parameter
subbox_per_rcHI = 0 # Neighborlist parameter.
######################
## CHEBYSHEV PARAMETERS
#######################
afmax = 0.35 # prefactor for FMAX
mfmax = 0.58 # exponent for FMAX
aract = 5.0 # prefactor for R_actual
mract = 0.63 # exponent for R_actual
Ncheb_incr = 2 #Increment in Number of Chebyshev terms (without changing lmax and lmin)
eigen_range = 5.0 #In percentage, used when increasing the width of lmax-lmin range
fd_err_range =100000 #For fd_err > fd_err_range, lmax-lmin range will increase by eigen_range (fd_err = fluctuation-dissipation error)
nrange_max = 500 #For nrange  > nrange_max, lmax-lmin range will increase by eigen_range. For example, if nrange_max is 1 then we provide just 2 opportunity to increment Ncheb, after that lmax-lmin range will increase by eigen_range %
Ncheb_evals_max = 100000 #This is the limit after which code will exit with a fatal error
fd_err_max = 0.01 #Error tolerance to satisfy FDT, 0.01 from Fixman's work
ns_interval = 50.0 # Time step interval at which we will check FDT error
ns_min_limit = 500 # Minimum limit of time step till we will check FDT error by default
#####################################################################################################
## SECTION 5: Simulation/Job length, time step etc...						    #
#####################################################################################################
## Relaxation time
#######################
model = 'Zimm' if (typeHI == 1) else 'Rouse' ## Use Rouse model for free draining and Zimm model for with HI (relaxation time)
## Rouse Model
tau_rouse = 0.5/(sin(0.5*pi/Nbpc)**2)##Eq 15.3-18 in Bird vol2
## Zimm Model
bzimm = 1.0 - (1.062233*(a**0.78)) ##Eq 15.4-28 in Bird vol2
sigmazimm = -0.895859*(a**0.78) ##Eq 15.4-29 in Bird vol2
azimm = 4.0*sin(0.5*pi/Nbpc)**2 *bzimm/(Nbpc**sigmazimm) ##Eq 15.4-27 in Bird vol2
tau_zimm = 2.0/azimm ##Eq 15.4-22 in Bird vol2
if model == 'Rouse':
    relaxation_time = tau_rouse
else:
    relaxation_time = tau_zimm
print 'relaxation time for a theta chain = ', relaxation_time
######################
## EQUILIBRATION
######################
equi_time = 3000.0
delT_Eq = 0.001                               # time step size
N_steps_Eq = int(equi_time/delT_Eq)            # Number of time steps
#N_steps_Eq = 10
data_interval_file_Eq = N_steps_Eq/10
data_interval_screen_Eq = N_steps_Eq/10     # After how many data points do you want to print on the screen?
######################
## PRODUCTION
######################
pro_time = 1000.0
delT_Pro = 0.001                             # time step size
N_steps_Pro = int(pro_time/delT_Pro)           # Number of time steps
#N_steps_Pro = 2
data_interval_file_Pro = N_steps_Pro/2000  # After how many data points do you want to print in data file?
#data_interval_file_Pro = N_steps_Pro/2
data_interval_screen_Pro = N_steps_Pro/10  # After how many data points do you want to print on the screen?
#data_interval_screen_Pro = N_steps_Pro/2
txtInterval = N_steps_Pro/10000
######################################################################################################
## SECTION 6: File names: I/O									     #
######################################################################################################
## SPRING TYPE
###################################
if springType == 0:
	spring = 'Hookean'
elif springType == 1:
        spring = 'FENE'
elif springType == 2:
        spring = 'ILC'
elif springType == 4:
        spring = 'FENE30'
elif springType == 5:
        spring = 'Fraenkel'  
else:
        spring = 'WLC'
###################################
## TRAJECTORY INFORMATION
###################################
if BATCH:
    import os
    run_number = int(os.environ['MY_VAR'])
    #run_number = rank + 1 
    print 'Trajectory#', run_number
else:
    run_number = 0
###################################
## FILE NAMES
###################################
#file_tag = "Gamma%sChi%sNb%sNc%sdt%s%sb%sphi%s" % (Gamma, Chi, Nbpc, Nc, delT_Pro, spring, b_para, phi)
file_tag = "Gamma%sChi%sNb%sNc%sdt%s%sb%sphi%sP%s" % (Gamma, Chi, Nbpc, Nc, delT_Pro, spring, b_para, phi, periodicity) # filetag for multisticker and telechelic polymers

## copying content of file1 to file2 corresponding to each run (only during production run)
file1 = "EQ%s.nc" % (file_tag)
file2 = "Run%i%s.nc" % (run_number, file_tag)
import shutil
shutil.copy2(file1, file2)

# INPUT
cdffile = '/short/g16/as6030/traj_f5cbycstr0.3/r%i%s.nc' % (run_number, file_tag) #The name of the .nc file
#cdffileEQ = "EQ%s.nc" % (file_tag)
cdffileEQ = "Run%i%s.nc" % (run_number, file_tag) # only during production run
if PREEQUILIBRATION:
    cdffile = cdffileEQ
###########################################END#######################################################
