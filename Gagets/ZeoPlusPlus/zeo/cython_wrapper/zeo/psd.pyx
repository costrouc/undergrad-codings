
from libcpp.string cimport string
from netstorage cimport AtomNetwork

def calcPoreSizeDistr(atmnet, highAcc_flag, probeChanRadius, probeRadius, 
        numSamples, excludePockets_flag, histFile, pointsFile, nodeRadiiFile,
        sphereDistFile, vis_flag, overlapCheck_flag):
    atmnet_copy = (<AtomNetwork?>atmnet).copy()
    c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    c_atmnetcp_ptr = (<AtomNetwork?>atmnet_copy).thisptr 
    cdef string chistFile = histFile
    cdef string cpntFile = pointsFile
    cdef string cndFile = nodeRadiiFile
    cdef string csphFile = sphereDistFile
    c_calcPoreSizeDistr (c_atmnetcp_ptr, c_atmnet_ptr, highAcc_flag,
              probeChanRadius,  probeRadius, numSamples, excludePockets_flag,
              chistFile, cpntFile, cndFile, csphFile, vis_flag, 
              overlapCheck_flag)

