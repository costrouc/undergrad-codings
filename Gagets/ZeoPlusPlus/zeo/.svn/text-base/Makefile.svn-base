# This Makefile compiles Zeo++ into an executable named network
# To see instructions of implementation please see README

#Use GNU C++ Compilier
CC = g++

#Flags to compile with
CFLAGS = -g

#Libraries to use while compiling
LIB = -lvoro++

# The relative of path of the main library source files
VOROINCLDIR=-I../../../Voro++/voro/trunk/src 
VOROLINKDIR=-L../../../Voro++/voro/trunk/src

#Object files to be created for network
OBJS=	main.o \
     network.o \
	networkinfo.o \
	networkio.o \
	instructions.o \
	networkanalysis.o \
	area_and_volume.o \
	segment.o \
	string_additions.o \
	symbcalc.o \
	graphstorage.o \
	networkstorage.o \
	voronoicell.o \
	channel.o \
	mindist.o \
	feature.o \
	geometry.o \
	holograms.o \
        grid.o \
	psd.o \
        sphere_approx.o \
        networkaccessibility.o \
        poreinfo.o

INCS= 	network.h \
	heap.h \
	general.h \
	networkstorage.h \
	networkinfo.h \
	mindist.h \
	networkio.h \
	geometry.h \
	graphstorage.h \
	graphstorage.h \
	instructions.h \
	voronoicell.h \
	channel.h \
	segment.h \
	networkanalysis.h \
	area_and_volume.h \
	holograms.h \
	string_additions.h \
	symbcalc.h \
	ray.h \
	feature.h \
	src/v_network.hh \
        grid.h \
        psd.h \
        sphere_approx.h \
        networkaccessibility.h \
        poreinfo.h

# List of the common source files

#Rich: added my file  Maciek added instuctions.cc and segment.cc

# List of executables
EXECUTABLES=network

# Makefile rules
all: network

network: ${OBJS} v_network.o ray.o
	@echo
	@echo Linking $@
	$(CC) $(VOROINCLDIR) $(VOROLINKDIR) -o $@ $(CFLAGS) ${OBJS} v_network.o ray.o $(LIB)

${OBJS}: ${INCS} $(@:%.o=%.cc) Makefile
	$(CC) $(VOROINCLDIR) -c $(@:%.o=%.cc) $(CFLAGS) 

ray.o: ray.cc ray.h
	$(CC) $(VOROINCLDIR) -c ray.cc $(CFLAGS) 

v_network.o: src/v_network.cc src/v_network.hh
	 $(CC) $(VOROINCLDIR) -c src/v_network.cc $(CFLAGS) 

clean:
	rm -f *.o *.err $(EXECUTABLES)

.PHONY: all clean
