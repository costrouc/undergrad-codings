# distutils: language = c++
# distutils: sources = ../networkstorage.cc
"""
Cython decloarations file for Zeo++ network storage section.
Declares Zeo++ ATOM_NETWORK, VORONOI_NETWORK classes and the associated 
python wrappers AtomNetwork, VoronoiNetwork.
"""

cdef extern from "../../networkstorage.h":
    cdef cppclass ATOM_NETWORK:
        ATOM_NETWORK() except +
        void copy(ATOM_NETWORK*)

    cdef cppclass VORONOI_NETWORK:
        VORONOI_NETWORK() except +
        VORONOI_NETWORK prune(double)

    cdef bint c_substituteAtoms "substituteAtoms"(ATOM_NETWORK*, ATOM_NETWORK*,
            bint, int*, bint)

    cdef bint c_fracSubstituteAtoms "fracSubstituteAtoms"(ATOM_NETWORK*, 
            ATOM_NETWORK*, bint, double, 
            int, int*, double*, bint)

cdef class AtomNetwork:
    """ 
    Cython wrapper declaration for Zeo++ ATOM_NETWORK class.
    Contains a pointer to ATOM_NETWORK and a flag denoting whether radisu
    for each atomic species is non-zero.
    """
    cdef ATOM_NETWORK* thisptr
    cdef bint rad_flag

cdef class VoronoiNetwork:
    """ 
    Cython wrapper declaration for Zeo++ VORONOI_NETWORK class.
    """
    cdef VORONOI_NETWORK* thisptr
