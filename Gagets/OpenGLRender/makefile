CC=gcc
CFLAGS=-c -Wall -g

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
LDFLAGS=
LIBDIR= -L/usr/lib64
LIB=-lX11 -lGL -lGLU -lglut -lm -lXmu
endif

ifeq ($(UNAME), Darwin)
LDFLAGS= -framework OpenGL -framework GLUT -framework Cocoa
LIBDIR=
LIB=
endif

DEFINED= -DGL_GLEXT_PROTOTYPES

SOURCES= meshobject.c \
	utils.c \
	phong.c \
	mouseinput.c \
	keyinput.c \
	shader.c \
	light.c \
	material.c \
	main.c \
	texture.c \
	lib/vec2.c \
	lib/vec3.c \
	lib/vec4.c \
	lib/mat2.c \
	lib/mat3.c \
	lib/mat4.c \
	lib/quat.c

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE= phong

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(LIBDIR) $(OBJECTS) $(LIB) -o $@

.c.o:
	$(CC) $(DEFINED) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) $(OBJECTS)

test: phong
	./phong plyfiles/bunny.ply
