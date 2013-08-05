cdef class VorCell:
    #cdef VOR_CELL *thiptr
    def __cinit__(self):
        self.thisptr = new VOR_CELL()

    def __dealloc__(self):
        del self.thisptr

cdef class BasicVCell:
    #cdef BASIC_VCELL *thisptr
    def __cinit__(self):
        self.thisptr = new BASIC_VCELL()

    def __dealloc__(self):
        del self.thisptr 
