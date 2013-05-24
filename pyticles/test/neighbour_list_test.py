import unittest

class NeighbourListTest(unittest.TestCase):
    """ 2. Test the neighbour lists. """
    """ The default neighbour list is a brute force list with no cutoff.
    """

    def test_nlist(self):
        """ Create a particle system with 3 particles and
            generate a brute force neighbour list.
            Check that the interparticle distances are correct.
        """
        import particles
        import neighbour_list
        n = 3
        p = particles.ParticleSystem(n,d=3,maxn=5)
        p.r[0,:] = (0.0,0.0,0.0) 
        p.r[1,:] = (1.0,0.0,0.0) 
        p.r[2,:] = (0.0,0.0,1.0) 
        nl = neighbour_list.NeighbourList(p)
        nl.build()
        nl.separations()
        k = nl.find_pair(0,1)
        self.assertEqual(nl.rij[k],1.0)

    def test_verlet_list(self):
        """ Create a 3 member particle system.
            Build, compress and compute separations for a verlet
            list.
            Ponder a rebuild - should be false.
            Move one particle a long way.
            Ponder a rebuild - should be true.

        """
        import particles
        import neighbour_list
        n = 3
        p = particles.ParticleSystem(n,d=3,maxn=5)
        p.r[0,:] = (0.0,0.0,0.0) 
        p.r[1,:] = (1.0,0.0,0.0) 
        p.r[2,:] = (0.0,0.0,1.0) 
        nl = neighbour_list.VerletList(p,cutoff=10,tolerance=2)
        nl.build()
        nl.compress()
        nl.separations()
        self.assertEqual(nl.rij[0],1.0)
        
        nl.ponder_rebuild()
        self.assertEqual(nl.rebuild_list,False)

        p.r[0,:] = (100.0,100.0,100.0)
        nl.ponder_rebuild()
        self.assertEqual(nl.rebuild_list,True)

    def test_fast_verlet(self):
        """ Create a 3 member particle system.
            Build, compress and compute separations for a verlet
            list.
            Ponder a rebuild - should be false.
            Move one particle a long way.
            Ponder a rebuild - should be true.

        """
        import particles
        import neighbour_list
        n = 3
        p = particles.ParticleSystem(n,d=3,maxn=5)
        p.r[0,:] = (0.0,0.0,0.0) 
        p.r[1,:] = (1.0,0.0,0.0) 
        p.r[2,:] = (0.0,0.0,1.0) 
        nl = neighbour_list.FastVerletList(p,cutoff=10,tolerance=2)
        nl.build()
        nl.compress()
        nl.separations()
        self.assertEqual(nl.rij[0],1.0)
        
        nl.ponder_rebuild()
        self.assertEqual(nl.rebuild_list,False)

        p.r[0,:] = (100.0,100.0,100.0)
        nl.ponder_rebuild()
        self.assertEqual(nl.rebuild_list,True)

 
if __name__=='__main__':
	unittest.main()

