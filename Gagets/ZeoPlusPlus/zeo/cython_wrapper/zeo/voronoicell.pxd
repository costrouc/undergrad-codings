cdef extern from "../../voronoicell.h":
    cdef cppclass VOR_CELL:
        VOR_CELL() except +

    cdef cppclass BASIC_VCELL:
        BASIC_VCELL() except +


cdef class VorCell:
    cdef VOR_CELL *thisptr

cdef class BasicVCell:
    cdef BASIC_VCELL *thisptr
