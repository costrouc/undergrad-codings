"""
Cython file defining methods for AtomNetwork and VoronoiNetowrk 
declared in netstorage.pxd file. 
"""

from libcpp.vector cimport vector

cimport zeo.netinfo, zeo.netio, zeo.network
from zeo.voronoicell cimport VorCell, BasicVCell, VOR_CELL, BASIC_VCELL
from zeo.network cimport visVoro

cdef class AtomNetwork:
    """
    Cython wrapper for Zeo++ ATOM_NETWORK class.
    Contains a pointer to ATOM_NETWORK and a flag denoting whether radisu
    for each atomic species is non-zero. The associated methods are 
    defined here.
    """
    #cdef ATOM_NETWORK* thisptr
    def __cinit__(self):
        self.thisptr = new ATOM_NETWORK()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    def copy(self):
        newatmnet = AtomNetwork()
        self.thisptr.copy(newatmnet.thisptr)
        return newatmnet

    @classmethod
    def read_from_CIF(cls, filename, rad_flag=True, rad_file=None):
        """
        Populates the ATOM_NETWORK from a CIF file.
        Calls Zeo++ readCIFFile function defined in networkio.cc.
        Arguments:
            filename: 
                Input CIF file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
        Returns:
            Instance of AtomNetwork
        """
        cdef char* c_rad_file = rad_file
        if rad_flag:
            if not rad_file:
                zeo.netinfo.zeo_initializeRadTable()
            else:       # rad_file is defined
                c_rad_file = rad_file
                zeo.netinfo.zeo_readRadTable(c_rad_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not zeo.netio.readCIFFile(c_filename, atmnet.thisptr, rad_flag):
            raise IOError
        atmnet.rad_flag = rad_flag
        return atmnet

    @classmethod
    def read_from_ARC(cls, filename, rad_flag=True, rad_file=None):
        """
        Populates the ATOM_NETWORK from a ARC file.
        Calls Zeo++ readARCFile function defined in networkio.cc.
        Arguments:
            filename: 
                Input ARC file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
        Returns:
            Instance of AtomNetwork
        """
        cdef char* c_rad_file = rad_file
        if rad_flag:
            if not rad_file:
                zeo.netinfo.zeo_initializeRadTable()
            else:       # rad_file is defined
                c_rad_file = rad_file
                zeo.netinfo.zeo_readRadTable(c_rad_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not zeo.netio.readARCFile(c_filename, atmnet.thisptr, rad_flag):
            raise IOError
        atmnet.rad_flag = rad_flag
        return atmnet

    @classmethod
    def read_from_CSSR(cls, filename, rad_flag=True, rad_file=None):
        """
        Populates the ATOM_NETWORK from a CSSR file.
        Calls Zeo++ readCSSRFile function defined in networkio.cc.
        Arguments:
            filename: 
                Input CSSR file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
        Returns:
            Instance of AtomNetwork
        """
        cdef char* c_rad_file
        if rad_flag:
            if not rad_file:
                zeo.netinfo.zeo_initializeRadTable()
            else:       # rad_file is defined
                c_rad_file = rad_file
                zeo.netinfo.zeo_readRadTable(c_rad_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not zeo.netio.readCSSRFile(c_filename, atmnet.thisptr, rad_flag):
            raise IOError
        atmnet.rad_flag = rad_flag
        return atmnet

    @classmethod
    def read_from_V1(cls, filename, rad_flag=True, rad_file=None):
        """
        Populates the ATOM_NETWORK from a V1 file.
        Calls Zeo++ readV1File function defined in networkio.cc.
        Arguments:
            filename: 
                Input V1 file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
        Returns:
            Instance of AtomNetwork
        """
        cdef char* c_rad_file = rad_file
        if rad_flag:
            if not rad_file:
                zeo.netinfo.zeo_initializeRadTable()
            else:       # rad_file is defined
                zeo.netinfo.zeo_readRadTable(c_rad_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not zeo.netio.readV1File(c_filename, atmnet.thisptr, rad_flag):
            raise IOError
        atmnet.rad_flag = rad_flag
        return atmnet

    def write_to_CSSR(self, filename):
        """
        Writes the ATOM_NETWORK to a CSSR file.
        Calls Zeo++ writeToCSSR function defined in networkio.cc.
        Arguments:
            filename: 
                Output CSSR file name.
        """
        cdef char* c_filename = filename
        if not zeo.netio.writeToCSSR(c_filename, self.thisptr):
            raise IOError

    def write_to_CIF(self, filename):
        """
        Writes the ATOM_NETWORK to a CIF file.
        Calls Zeo++ writeToCIF function defined in networkio.cc.
        Arguments:
            filename: 
                Output CIF file name.
        """
        cdef char* c_filename = filename
        if not zeo.netio.writeToCIF(c_filename, self.thisptr):
            raise IOError

    def write_to_V1(self, filename):
        """
        Writes the ATOM_NETWORK to a V1 file.
        Calls Zeo++ writeToV1 function defined in networkio.cc.
        Arguments:
            filename: 
                Output V1 file name.
        """
        cdef char* c_filename = filename
        if not zeo.netio.writeToV1(c_filename, self.thisptr):
            raise IOError

    def write_to_XYZ(self, filename, supercell_flag, 
                     is_duplicate_perimeter_atoms):
        """
        Writes the ATOM_NETWORK to an XYZ file.
        Calls Zeo++ writeToXYZ function defined in networkio.cc.
        Arguments:
            filename: 
                Output XYZ file name.
            supercell_flag:
                Flag denoting whether to write 2x2x2 supercell.
            is_duplicate_perimeter_atoms:
                Flag denoting whether perimeter atoms need to be replicated.
        """
        cdef char* c_filename = filename
        if not zeo.netio.writeToXYZ(c_filename, self.thisptr, supercell_flag, 
                is_duplicate_perimeter_atoms):
            raise IOError

    def write_to_VTK(self, filename):
        """
        Writes the boundary of unit cell within the ATOM_NETWORK to a VTK file.
        Calls Zeo++ writeToVTK function defined in networkio.cc.
        Arguments:
            filename: 
                Output VTK file name.
        """
        cdef char* c_filename = filename
        if not zeo.netio.writeToVTK(c_filename, self.thisptr):
            raise IOError

    def write_to_MOPAC(self, filename, supercell_flag):
        cdef char* c_filename = filename
        if not zeo.netio.writeToMOPAC(c_filename, self.thisptr, supercell_flag):
             raise IOError

    def calculate_free_sphere_parameters(self, filename):
        vornet = self.perform_voronoi_decomposition()
        cdef char* c_fname = filename
        vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
        zeo.network.calculateFreeSphereParameters(vornet_ptr, c_fname, False)

    def perform_voronoi_decomposition(self, saveVorCells=False):
        """
        Performs voronoi decomposition for the ATOM_NETWORK to analyze
        void space and generate voronoi nodes, edges and faces
        Calls Zeo++ performVoronoiDecomp function defined in network.cc.
        Arguments:
            saveVorCells (optional): 
                Flag to denote whether to save the VorCells.
                Reserved for future use, so ignore this.
        Returns:
            Instance of VoronoiNetwork
        """
        vornet = VoronoiNetwork()  
        cdef vector[VOR_CELL] vcells
        cdef vector[BASIC_VCELL] bvcells
        print self.rad_flag
        if not zeo.network.performVoronoiDecomp(self.rad_flag, self.thisptr, 
                vornet.thisptr, &vcells, saveVorCells, &bvcells):
            raise ValueError # Change it to appropriate error
        cdef int N
        vorcelllist = []
        bvcelllist = []
        # define copy methods for BASIC_VCELL and VOR_CELL methods
        for i in range(vcells.size()):
            pass
            #vorcell = VorCell()
            #vorcell.thisptr = &(vcells[i])
            #vorcelllist.append(vcells[i])
        for i in range(bvcells.size()):
            pass
            #basicvcell = BasicVCell()
            #basicvcell.thisptr = &(bvcells[i])
            #bvcelllist.append(bvcells[i])
        return vornet


cdef class VoronoiNetwork:
    #cdef VORONOI_NETWORK* thisptr
    def __cinit__(self):
        self.thisptr = new VORONOI_NETWORK()
    def __init__(self):
        pass
    def __dealloc__(self):
        del self.thisptr
    def prune(self, minRad):
        cdef VORONOI_NETWORK newcvornet = self.thisptr.prune(minRad)
        newvornet = VoronoiNetwork()
        newvornet.thisptr = &newcvornet
        return newvornet

    def analyze_writeto_XYZ(self, name, double probeRad, atmnet, 
            int shift_x=0, int shift_y=0, int shift_z=0):
        cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
        cdef char* cname = name
        visVoro(name, probeRad, shift_x, shift_y, shift_z, self.thisptr, 
                c_atmnetptr)


def substitute_atoms(atmnet, substituteFlag, radialFlag):
    cdef int substitutionNo[1]
    atmnet_copy = AtomNetwork()
    c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    if not c_substituteAtoms(c_atmnet_ptr, atmnet_copy.thisptr, substituteFlag,
            substitutionNo, radialFlag):
        raise ValueError
    subNo = substitutionNo[0]
    return atmnet_copy, subNo



