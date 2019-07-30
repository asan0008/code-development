### Time-stamp: <PolymerModel.py 19:21, 09 Dec 2009 by P. Sunthar>

"""
Coarse grained polymer model objects.

This module contains definitions for coarse grained representation of polymers.
Models such as Bead Spring, Bead Rod (not implemented).

"""

import cluster_customization

__docformat__ = 'epytext'

from MMTK.ChemicalObjects  import *
from MMTK import Bonds
from Scientific.Geometry import Vector
from MMTK.ForceFields.ForceField import ForceField
#execfile('Input.py')
from Input import *
## Alternate Random number generator that accepts seed  ###
## seed can be set in calling program by
## from numpy.random import *
## seed(1234)

## MMTK.Random always uses random seed, not helpful for debugging
## from MMTK.Random import *

import numpy.random 
def gaussian(mean, std, shape=None):
    return numpy.random.normal(mean,std,shape)

def randomPointInBox(a, b = None, c = None):
    """
    @param a: the edge length of a box along the x axis
    @type a: C{float}
    @param b: the edge length of a box along the y axis (default: a)
    @type b: C{float}
    @param c: the edge length of a box along the z axis (default: a)
    @type c: C{float}
    @returns: a vector drawn from a uniform distribution within a
              rectangular box with edge lengths a, b, c.
    @rtype: C{Scientific.Geometry.Vector}
    """
    if b is None: b = a
    if c is None: c = a
    x = numpy.random.uniform(-0.5*a, 0.5*a)
    y = numpy.random.uniform(-0.5*b, 0.5*b)
    z = numpy.random.uniform(-0.5*c, 0.5*c)
    return Vector(x, y, z)

def centerPointInBox(a, b = None, c = None):
    """
    @param a: the edge length of a box along the x axis
    @type a: C{float}
    @param b: the edge length of a box along the y axis (default: a)
    @type b: C{float}
    @param c: the edge length of a box along the z axis (default: a)
    @type c: C{float}
    @returns: a vector drawn from a uniform distribution within a
              rectangular box with edge lengths a, b, c.
    @rtype: C{Scientific.Geometry.Vector}
    """
    if b is None: b = a
    if c is None: c = a
    x = 0.0
    y = 0.0
    z = 0.0
    return Vector(x, y, z)


class BeadCluster(AtomCluster):
    def __init__(self, NBeads=2, name='chain', element='c', atomdict = {}):
        if NBeads < 2:
            raise ValueError('Number of Beads must be >=2')

        beadlist = [Atom(element, name='B' + str(x), **atomdict)
                    for x in range(NBeads)]
        
        AtomCluster.__init__(self, beadlist, name=name)
        self.NBeads = len(self.atoms)  

    def clearBondAttributes(self):
        for a in self.atoms:
            a.clearBondAttribute()
        
    def bondedTo(self, atom):
        return self.bonds.bondedTo(atom)

from MMTK.Bonds import Bond

## This should tally with the numbers used in
## MMTK_entropic_spring.c: entropic_spring_evaluator()

SpringTypes = { 'Hookean' : 0,
                'FENE'    : 1,
                'ILC'     : 2,
                'WLC'     : 3,
                'FENE30'  : 4,
                'Fraenkel': 5
                }
    

class Spring(Bond):
    def __init__(self, (a1, a2), sptype, b_parameter):
        """ (a1, a2) : Atoms that will be joined by the spring
            sptype   : string that has to be convertable by SpringTypes dictionary
            b_parameter: FIXME
        """
        Bonds.Bond.__init__(self, (a1, a2))
        try:
            self.sptype = SpringTypes[sptype]
        except KeyError:
            raise ValueError( sptype +  'spring force not defined.')
        
        if not (isinstance(b_parameter, float)):
            raise ValueError("b_parameter must be a floating point number.")
        self.b_parameter = b_parameter

from MMTK import Utility
from numpy import *

class LinearChain(BeadCluster):
    
    def __init__(self, NBeads=2, sptype='FENE',
                 b_parameter=2, name='chain', element='c', atomdict = {}):
        """Constructs a LinearChain of NBeads that are joined by springs of type sptype."""

        BeadCluster.__init__(self, NBeads, name, element, atomdict)
        self.bonds = Bonds.BondList([])
         
        for b in range(NBeads-1):
            self.bonds.append(Spring((self.atoms[b], self.atoms[b+1]), sptype, b_parameter))

        self.bonded_options = {'bonds': True, 'bond_angles': False,
                               'dihedrals': False, 'impropers': False}

        self.initialiseBeadPositions(b_parameter,sptype)

    def initialiseBeadPositions(self,b_parameter,sptype):
        if sptype == 'Hookean':
            rvec = cumsum(gaussian(0,1, (self.NBeads-1, 3)),axis=0)
        else:
            # for all other spring forces treat them to lie in a box
            # centered at origin. This means the max distance is
            # diagonal/2 = sqrt(b) giving a = 2 sqrt(b/3). This is a
            # quick approximation to avoid overstreched springs.
            a = 2*sqrt(b_parameter/3)
            #a = sqrt(b_parameter/3) # reduced the maximum distance between two adjacent beads to avoid overstretching
            rvec = cumsum(array([randomPointInBox(a).array
                                  for x in range(self.NBeads-1)]), axis=0)
            
        b = self.bonds[0]
        b.a1.translateTo(Vector(0,0,0))
                
        for (b, r) in zip(self.bonds, rvec):
            b.a2.translateTo(Vector(r))


from MMTK.ForceFields.BondedInteractions import BondedForceField

class EntropicSpringForceField(BondedForceField):
    """Adapted from ForceFields.MMForceField.MMBondedForceField."""
    def __init__(self):
        BondedForceField.__init__(self, "Entropic Spring Force")
        self.arguments =  () 
        # used in Trajectory, to store universe.description() in skeleton
        # must be a tuple.

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        return BondedForceField.evaluatorParameters(self, universe,
                                                    subset1, subset2,
                                                    global_data)

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        from MMTK_entropic_spring import EntropicSpringTerm
        param = self.evaluatorParameters(universe, subset1, subset2,
                                         global_data)
        bonds = param['harmonic_distance_term']
        if bonds:
            indices = N.array(map(lambda b: b[:2], bonds))
            sptypes = N.array(map(lambda b: b[2], bonds))
            b_params = N.array(map(lambda b: b[3], bonds))
            en_term = EntropicSpringTerm(universe._spec,
                                           indices, sptypes, b_params)
        return [en_term]

    ## Functions required to be defined in BondedInteractions.BondedForceField
    ## adapted from MMForceField.py
    
    def addBondTerm(self, data, bond, object, global_data):
        a1 = bond.a1
        a2 = bond.a2
        i1 = a1.index
        i2 = a2.index

        ptyp = bond.sptype
        pval = bond.b_parameter

        data.add('bonds', (i1, i2, ptyp, pval)) #bond indices spring param passed to c function


    def bonds(self, global_data):
        return self.data.get('bonds')

    ## Other functions that are required or recommended should raise
    ## the same AttributeError here (as the base class) as they are
    ## undefined for this force field

from MMTK import Dynamics, Features, Trajectory, ParticleProperties

# adds this MMTK_brownian1 to the searchpath
import sys
sys.path.append(cluster_customization.MMTK_brownian_path)

if cluster_customization.name == "sungrid":
    sys.path.append("/nfs/monash/home/jain/usr/lib/python")
    sys.path.append("/opt/sw/atlas-3.8.1/lib")

from MMTK_brownian1 import integrateBD
class BrownianIntegrator(Dynamics.Integrator):
    """ Brownian Dynamics Integrator """
    def __init__(self, universe, **options):
        Dynamics.Integrator.__init__(self, universe, options)
        # Supported features: none for the moment, to keep it simple
        self.features = []

    def __call__(self, **options):
        # Process the keyword arguments
        self.setCallOptions(options)
        # Check if the universe has features not supported by the integrator
        Features.checkFeatures(self, self.universe)
        # Get the universe variables needed by the integrator
        configuration = self.universe.configuration()

        # Get the radii of beads. First check if a keyword argument
        # 'radius' was given to the integrator (radius=<value>). In
        # this case, its value can be a ParticleScalar or a plain
        # number (used for all atoms). If no such argument is given,
        # collect the values of the attribute 'radius' from all atoms
        # (default is zero).
        try:
            radii = self.getOption('radius')
        except ValueError:
            radii = self.universe.getParticleScalar('radius')
        if not ParticleProperties.isParticleProperty(radii):
            var = ParticleProperties.ParticleScalar(self.universe)
            var.array[:] = radii
            radii = var
        if max(radii) != min(radii):
            raise ValueError('Multiple radii not implemented.')

        # for ES
        charges = array([getattr(a, chargename) if hasattr(a, chargename) else 0 for a in self.universe.atomList()])

        #lx = self.universe.boxToRealCoordinates(Vector(1., 0., 0.)).length()
        #ly = self.universe.boxToRealCoordinates(Vector(0., 1., 0.)).length()
        #lz = self.universe.boxToRealCoordinates(Vector(0., 0., 1.)).length()

        # Construct a C evaluator object for the force field, using
        # the specified number of threads or the default value
        nt = self.getOption('threads')
        evaluator = self.universe.energyEvaluator(threads=nt).CEvaluator()
        print 'You are in Polymer Model' 
        # Run the C integrator
        integrateBD(self.universe,
                    configuration.array, 
                    radii.array,
                    charges,
                    evaluator,
                    self.getOption('delta_t'),
                    self.getOption('start_time'),
                    self.getOption('tStrain'),
                    self.getOption('iStrain'),
                    self.getOption('steps'),
                    self.getOption('fixL1x'),
                    self.getOption('fixL1y'),
                    self.getOption('fixL2x'),
                    self.getOption('fixL2y'),
                    cbycstar,
                    typeHI,
                    Nbpc,
		    springType,
		    b_para,
		    SHEAR,
		    ELONGATION,
		    EQinEQ,
		    gam_dot,
                    eps_dot,
		    taup_MF,
                    taup_SH,
		    theta_MF,
                    M,
                    factor1,
                    factor2,
                    subbox_per_rcHI,
                    afmax,
                    mfmax,
                    aract,
                    mract,
                    Ncheb_incr,
                    eigen_range,
                    fd_err_range,
                    nrange_max,
                    Ncheb_evals_max,
                    fd_err_max,   
                    ns_interval,
                    ns_min_limit,
                    self.getOption('EQ'),
                    self.getOption('weiner_seed1'),
                    self.getOption('weiner_seed2'),
                    run_number,
                    self.getActions(),
                    typeEV,
                    subbox_per_rcEV,
                    EV_cutoff,
                    lj1,
                    lj2,
                    dst,
                    zst,
                    typeAS,
                    periodicity,
                    subbox_per_rcAS,
                    AS_cutoff1,
                    AS_cutoff2,  
                    Lja1,
                    Lja2,
                    phi,
                    alphaAS,
                    betaAS, 
                    typeDH,
                    subbox_per_rcDH,
                    DH_cutoff,
                    Lbstar,
                    Ldstar,
                    typeES,
                    subbox_per_rcES,
                    energy_factorES,
                    self.getOption('intval'),
		    txtInterval,
		    RA,
                    timeStepFactor,
                    N_wait_rampup,
                    Q_para  
                    )


