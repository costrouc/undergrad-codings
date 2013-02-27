CC=g++
CFLAGS=-c -Wall -g
LDFLAGS=
SOURCES=Octree.cpp \
	OctreeNode.cpp \
	Triple.cpp \
	Point3d.cpp \
	Vector3d.cpp \
	Ray.cpp \
	Sphere.cpp \
	Plane.cpp \
	main.cpp

OBJECTS=$(SOURCES:.cpp=.o)
INCLUDES=$(SOURCES:.cpp=.h)

EXECUTABLE=octree

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o: 
	$(CC) $(CFLAGS) $< -o $@

clean:
	@echo Removing object files...
	rm $(OBJECTS)