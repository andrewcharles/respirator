import unittest

class KernelTest(unittest.TestCase):
    """ 1. Test the smooth particle kernel functions. """

    def test_c_kernel(self):
        import c_test
        c_test.kernel_test()

    def test_py_kernel(self):
        import spkernel
        w, dwdx = spkernel.lucy_kernel(0.0,[0.0,0.0,0.0],2.0)
        print 'Python kernel at zero distance ',w

if __name__=='__main__':
	unittest.main()
