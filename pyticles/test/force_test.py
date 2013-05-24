import unittest

class ForceTest(unittest.TestCase):
    """ Test the sph force subroutine. """

    def test_spam_complete_force(self):
        """ """
        pass

    def test_spam_force(self):
        import particles
        import neighbour_list
        import properties
        import forces
        n = 2
        p = particles.SmoothParticleSystem(n,d=3,maxn=5)
        p.r[0,:] = (0.0,0.0,0.0) 
        p.r[1,:] = (1.0,0.0,0.0) 
        nl = neighbour_list.VerletList(p,cutoff=10,tolerance=2)
        nl.build()
        nl.compress()
        nl.separations()
        properties.spam_properties(p,nl)
        f = forces.SpamForce(p,nl)
        f.apply()

    def test_fortran_collision(self):
        import particles
        import neighbour_list
        import properties
        import forces
        n = 2
        p = particles.SmoothParticleSystem(n,d=3,maxn=5)
        p.r[0,:] = (0.0,0.0,0.0) 
        p.r[1,:] = (1.0,0.0,0.0) 
        p.v[0,:] = (1.0,0.0,0.0) 
        p.v[1,:] = (0.0,0.0,0.0) 
        nl = neighbour_list.VerletList(p,cutoff=10,tolerance=2)
        nl.build()
        nl.compress()
        nl.separations()
        f = forces.FortranCollisionForce(p,nl,cutoff=2.0)
        f.apply()
        print p.v[0,:]
        print p.v[1,:]

    def test_python_collision(self):
        import particles
        import neighbour_list
        import properties
        import forces
        n = 2
        p = particles.SmoothParticleSystem(n,d=3,maxn=5)
        p.r[0,:] = (0.0,0.0,0.0) 
        p.r[1,:] = (1.0,0.0,0.0) 
        p.v[0,:] = (1.0,0.0,0.0) 
        p.v[1,:] = (0.0,0.0,0.0) 
        nl = neighbour_list.VerletList(p,cutoff=10,tolerance=2)
        nl.build()
        nl.compress()
        nl.separations()
        f = forces.CollisionForce3d(p,nl,cutoff=2.0)
        f.apply()
        print p.v[0,:]
        print p.v[1,:]



if __name__=='__main__':
	unittest.main()
