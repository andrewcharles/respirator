cimport c_properties

def kernel_test():
    c_kernel_test()

cdef void c_kernel_test():
    cdef float r, dx[3], h, w, dwdx[3]
    r = 0.0
    dx[0] = 0.0
    dx[1] = 0.0
    dx[2] = 0.0
    h = 2.0
    w = 0.0
    dwdx[0] = 0.0
    dwdx[1] = 0.0
    dwdx[2] = 0.0
    c_properties.lucy_kernel_3d(r,dx,h,&w,dwdx)
    print 'C kernel at zero distance',w

