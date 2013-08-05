# distutils: language = c++
# distutils: sources = ../graphstorage.cc

from netstorage cimport VORONOI_NETWORK

cdef extern from "../../graphstorage.h":
    cdef cppclass DIJKSTRA_NETWORK:
        DIJKSTRA_NETWORK() except +
cdef extern from "../../graphstorage.h" namespace "DIJKSTRA_NETWORK":
    cdef void buildDijkstraNetwork(VORONOI_NETWORK*, DIJKSTRA_NETWORK*)


cdef class DijkstraNetwork:
    cdef DIJKSTRA_NETWORK* thisptr
