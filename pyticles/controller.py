""" Objects for mutating particle systems. Pair potentials, body forces
    and even boundary conditions can all be modelled as controllers.

"""

import numpy

class Controller:
    """ Mutates one ParticleSystem.
    """

    def __init__(self,p):
        """ Initialise the controller with a reference to the
            particle system it is to mutate.

            p -- ParticleSystem object to be mutated.

        """
        self.p = p
        return
    
    def apply(self):
        """All controllers must have an apply method."""
        print "Do nothing"
        return

class CouplingController:
    """ Mutates one or more ParticleSystems.
    """

    def __init__(self):
        """ Initialise the controller with a reference to the
            particle system it is to mutate.


        """
        self.groups = []
        return

    def bind_particles(self,p):
        """
            p -- ParticleSystem object to be mutated.
        """
        self.groups.append(p)

    def apply(self):
        """All controllers must have an apply method."""
        print "Do nothing"
        return


#class PairForce(Force):
#    """ A force that operates between pairs of particles.
#        This one needs a data structure to keep track of the
#        pairs.
#    """
#
#    def __init__(self,nl):
#        """ We only need to initialise with the neighbour list
#            because the neighbour list refers to the particles
#            it points to. There may be a clever way to implement
#            this so that the pair force routines do not need to
#            know if the two particles are from different systems
#            or not.
#        """


