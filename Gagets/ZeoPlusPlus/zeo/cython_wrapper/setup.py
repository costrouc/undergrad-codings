#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
#from Cython.Build import cythonize

includedirs=["../../../../Voro++/voro/trunk/src"]
libdirs = ["../../../../Voro++/voro/trunk/src"]
netstorage_srcfiles = [
        'zeo/netstorage.pxd', 'zeo/netstorage.pyx', '../networkstorage.cc', 
        '../mindist.cc', '../geometry.cc', '../networkinfo.cc',
        '../networkio.cc', '../grid.cc', '../symbcalc.cc',
        '../string_additions.cc', 'zeo/voronoicell.pyx', '../voronoicell.cc', 
        '../networkanalysis.cc', '../graphstorage.cc', '../area_and_volume.cc',
        'zeo/network.pxd', '../network.cc', '../src/v_network.cc',
        '../networkaccessibility.cc', '../channel.cc'
        ]
netinfo_srcfiles = ['zeo/netinfo.pyx', '../networkinfo.cc']
netio_srcfiles = [
        'zeo/netio.pyx', '../networkio.cc', 'zeo/netinfo.pyx', 
        '../networkinfo.cc', 'zeo/string_add.pxd', '../string_additions.cc', 
        '../grid.cc', '../mindist.cc', '../symbcalc.cc', 
        '../networkstorage.cc', '../geometry.cc'
        ]
graphstorage_srcfiles = ['zeo/graphstorage.pyx', '../graphstorage.cc']
psd_srcfiles = ['zeo/psd.pyx', '../psd.cc']
voronoicell_srcfiles = [
        'zeo/voronoicell.pyx', '../voronoicell.cc'
        ]
channel_srcfiles = ['zeo/channel.pyx', '../channel.cc']
highaccuracy_srcfiles = [
        'zeo/high_accuracy.pyx', '../sphere_approx.cc', '../networkstorage.cc',
        '../networkinfo.cc', '../mindist.cc', '../geometry.cc'
        ]
areavol_srcfiles = [
        'zeo/area_volume.pyx', '../area_and_volume.cc'#,
        #'../networkinfo.cc', '../mindist.cc', '../geometry.cc'
        ]
setup(
    name = 'zeo',
    version = '0.1',
    description = "Python interface to Zeo++",
    url = "",
    author = "Maciej Haranzcyk",
    license = "",
    cmdclass = {'build_ext':build_ext},
    ext_modules = [Extension("zeo.voronoicell", voronoicell_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'], 
                             library_dirs = libdirs,
                             language = 'c++'),
                   Extension("zeo.netstorage", netstorage_srcfiles, 
                             include_dirs=includedirs,
                             libraries = ['voro++'], 
                             library_dirs = libdirs,
                             language = 'c++'),
                   Extension("zeo.netinfo", netinfo_srcfiles,
                             language = 'c++'),
                   Extension("zeo.netio", netio_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             language = 'c++'),
                   Extension("zeo.graphstorage", graphstorage_srcfiles,
                             language = 'c++'),
                   Extension("zeo.psd", psd_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             language = 'c++'),
                   Extension("zeo.channel", channel_srcfiles,
                             language = 'c++'),
                   Extension("zeo.high_accuracy", highaccuracy_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             language = 'c++'),
                   Extension("zeo.area_volume", areavol_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             language = 'c++'),
                             ]
)
