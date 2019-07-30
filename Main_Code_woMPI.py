execfile('Input.py')
from MMTK import Vector
###################################################
# INITIAL BASIS VECTORS OF THE SIMULATION BOX
# MIXED FLOW
B10x = 1.0
B10y = 0.0
B10z = 0.0
B20x = 0.0
B20y = 1.0
B20z = 0.0
B30x = 0.0
B30y = 0.0
B30z = 1.0

if(SHEAR==1 and ELONGATION==1):
	B10x = cos(theta_MF) - ((0.5*gam_dot/eps_dot)*sin(theta_MF))
	B10y = sin(theta_MF)
	B10z = 0.0
	B20x = -sin(theta_MF) - ((0.5*gam_dot/eps_dot)*cos(theta_MF))
	B20y = cos(theta_MF)
	B20z = 0.0
	B30x = 0.0
	B30y = 0.0
	B30z = 1.0
# ONLY ELONGATIONAL FLOW
if(ELONGATION==1 and SHEAR==0):
	B10x = cos(theta_MF)
	B10y = sin(theta_MF)
	B10z = 0.0
	B20x = -sin(theta_MF)
	B20y = cos(theta_MF)
	B20z = 0.0
	B30x = 0.0
	B30y = 0.0
	B30z = 1.0
# ONLY SHEAR FLOW TODO: Check the following equations for "only shear flow" case
if(SHEAR==1 and ELONGATION==0):
	B10x = 1.0
	B10y = 0.0
	B10z = 0.0
        # For -45 to +45 system
        B20x = -1.0
	# For 0 to +45 system 
	#B20x = 0.0
	B20y = 1.0
	B20z = 0.0
	B30x = 0.0
	B30y = 0.0
	B30z = 1.0
if(EQinEQ == 1):
        B10x = 1.0
        B10y = 0.0
        B10z = 0.0
        B20x = 0.0
        B20y = 1.0
        B20z = 0.0
        B30x = 0.0
        B30y = 0.0
        B30z = 1.0
##################################################
# INITIAL BOX SIZES
L10 = Vector(Lx*B10x, Lx*B10y, Lx*B10z)
L20 = Vector(Ly*B20x, Ly*B20y, Ly*B20z)
L30 = Vector(Lz*B30x, Lz*B30y, Lz*B30z)
print 'L10', L10
print 'L20', L20
print 'L30', L30
vector2 = [L20[1]*L30[2] - L20[2]*L30[1],
         L20[2]*L30[0] - L20[0]*L30[2],
         L20[0]*L30[1] - L20[1]*L30[0]]
volumeBox = sum(p*q for p,q in zip(L10, vector2))
print "Volume1 = ", volumeBox
print "Volume2 = ", Lx*Ly*Lz
#phi_MF = atan(L20[1]/L20[0])
#theta1_MF = atan(L10[1]/L10[0])
#print 'theta1_MF = ', theta1_MF
if(SHEAR==1 and ELONGATION==0):
	taup_SH = abs(2.0*L20[0]/(gam_dot*L20[1]))
else:
	taup_SH = 0.0 #Dummy value for taup_SH (# Irrelevant for extensional or mixed flow)
#Storing "fixed" box dimensions
fixL1x = L10[0]
fixL1y = L10[1]
fixL2x = L20[0]
fixL2y = L20[1]
print "fixL1x = ", fixL1x
print "fixL1y = ", fixL1y
print "fixL2x = ", fixL2x
print "fixL2y = ", fixL2y
###################################################
print '############################################'
print 'INPUT PARAMETERS'
print '############################################'
print 'c/c* = ', cbycstar
print 'Nbpc = ', Nbpc
print 'bead radius = ', a
print 'Rg_Rouse = ', Rg
print 'Re_approx = ', Re
print 'Lx(Orthorhombic) = ', Lx
print 'Ly(Orthorhombic) = ', Ly
print 'Lz(Orthorhombic) = ', Lz
print 'Lx(Parallelepipedic) = ', L10.length()
print 'Ly(Parallelepipedic) = ', L20.length()
print 'Lz(Parallelepipedic) = ', L30.length()
print 'Nc = ', Nc
if springType != 1:
	print 'b = ', b_para
if(SHEAR==1 and ELONGATION==0): #ONLY SHEAR
	print 'Gamma dot = ',gam_dot
	print 'Lattice strain period_SH = ', taup_SH
if(SHEAR==0 and ELONGATION==1): #ONLY ELONGATION
	print 'Epsilon dot = ',eps_dot
        print 'Lambda_p_MF = ', lamdap_MF
        print 'N11 = ', N11
        print 'N12 = ', N12
	print 'Henky strain_MF = ', eps_p_MF
	print 'Lattice strain period_MF = ', taup_MF
	print 'Magic angle in rad (MF) = ', theta_MF
if(SHEAR==1 and ELONGATION==1): #MIXED FLOW
	print 'Gamma dot = ',gam_dot
	print 'Epsilon dot = ',eps_dot
	print 'Lambda_p_MF = ', lamdap_MF
	print 'Henky strain_MF = ', eps_p_MF
	print 'Lattice strain period_MF = ', taup_MF
	print 'Magic angle in rad (MF) = ', theta_MF
print 'tau_rouse = ', tau_rouse
print 'tau_zimm = ', tau_zimm
print 'typeHI = ', typeHI
print 'typeEV = ', typeEV
print 'typeAS = ', typeAS
print 'ASPeriodicity = ', periodicity
print 'AS_cutoff1= ', AS_cutoff1
print 'AS_cutoff2= ', AS_cutoff2
if (typeEV == 2): # Narrow Gaussian
	print 'EV_cutoff = ', EV_cutoff
	print 'd* = ', dst
	print 'z* = ', zst
        print 'EVcutoff set originally = ', EV_cutoff
print 'equi_time = ', equi_time
print 'delT_Eq = ', delT_Eq
print 'N_steps_Eq = ', N_steps_Eq
print 'pro_time = ', pro_time
print 'delT_Pro = ', delT_Pro
print 'N_steps_Pro = ', N_steps_Pro
print '############################################'

### Use the same seed for debugging
from numpy.random import *
import random
#seed(1234) # seed for placing beads randomly in a chain
append_filename = 'append_link_%s.dat' % run_number
#seed_generation = 235485
seed_generation = int(random.randint(10000,999999))
seed(seed_generation) 
print 'seed_generation = ', seed_generation

execfile('sysCreate.py')
from MMTK import *
from MMTK.Trajectory import Trajectory, TrajectoryOutput, StandardLogOutput
yzGap = L30.length()/Nc
xGap = 0.9*sqrt(b_para)
#trajectory = Trajectory(universe, cdffile, "a", "Full_Configurations", double_precision = True)
#full = Trajectory(universe, cdffileEQ)  # refer MMTK library
for i in range(1):
    if PREEQUILIBRATION:
        ## Placing RCOM of each chain randomly between -L/2  to +L/2
        for m  in universe.objectList():
	    rcom = randomPointInBox(L10[0], L20[1], L30.length())#Uniform distribution
            m.translateTo(rcom)
    else:
        full = Trajectory(universe, cdffileEQ)  # refer MMTK library
        universe.setFromTrajectory(full, len(full)-1) # refer MMTK
        print '####################'
    print universe.configuration().array
    print 'run no. = ', i
    weiner_seed1 = int(random.randint(1,2**24 - 3))
    weiner_seed2 = int(random.randint(1,2**24 - 3))
    #weiner_seed1 = 34654656
    #weiner_seed2 = 56657675
    if weiner_seed1 % 2 == 0:
            weiner_seed1 += 1
    initialtime = 0.0
    timeStrain = 0.0
    idStrain = 1

    print 'weiner_seed1 = ', weiner_seed1
    print 'weiner_seed2 = ', weiner_seed2

    universe.setForceField(EntropicSpringForceField())
    trajectory = Trajectory(universe, cdffile, "w", "Full_Configurations", double_precision = True)
    # Create integrator
    integrator = BrownianIntegrator(universe) # refer MMTK
    print 'You are in Main_Code.py'

    # Equilibration
    integrator(steps=N_steps_Eq, start_time = initialtime, tStrain = timeStrain, iStrain = idStrain, delta_t=delT_Eq, fixL1x = fixL1x, fixL1y = fixL1y, fixL2x = fixL2x, fixL2y = fixL2y, EQ = 1, weiner_seed1 = weiner_seed1, weiner_seed2 = weiner_seed2, intval=data_interval_file_Eq,
                    actions=[StandardLogOutput(data_interval_screen_Eq)])

    # get random seed from equilibration run
    file = open(append_filename, "r")
    weiner_seed1 = float(file.readline())
    weiner_seed2 = float(file.readline())
    dummyTime = float(file.readline())
    print 'weiner_seed1 = ', weiner_seed1
    print 'weiner_seed2 = ', weiner_seed2
    #timeStrain = float(file.readline())
    #idStrain = int(file.readline())
    timeStrain = 0
    idStrain = 1
    # Production
    integrator(steps=N_steps_Pro, start_time = initialtime, tStrain = timeStrain, iStrain = idStrain, delta_t=delT_Pro, fixL1x = fixL1x, fixL1y = fixL1y, fixL2x = fixL2x, fixL2y = fixL2y, EQ = 0, weiner_seed1 = weiner_seed1, weiner_seed2 = weiner_seed2, intval=data_interval_file_Pro,
                    actions = [ # Write every ## step to the trajectory file.
                               TrajectoryOutput(trajectory, ("time", "energy",
                                                             "configuration","gradients" ),
                                                0, None, data_interval_file_Pro),
                               # Log output to screen every ## steps.
                               StandardLogOutput(data_interval_screen_Pro)])
trajectory.close()
'''
#POSTPROCESSING
if SERIES==False or NJOBS==NJOB:
    cdffile=cdffilePRO
    import postprocessing
    postprocessing.process(cdffile, file_output, All_Prop, SF)
'''
