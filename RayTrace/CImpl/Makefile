# Makefile for ray 1
CC = gcc
CFLAGS = -g -Wall

OBJS= ray1.o \
      object.o \
      shape.o \
      world.o \
      window.o \
      plane.o \
      sphere.o \
      image.o \
      raytrace.o   \
      list.o \
      getColor.o \
      triple.o \
      polygon.o \
      light.o \
      illuminate.o \
      rectangle.o \
      triangle.o

INCS= ray.h \
      image.h \
      light.h \
      list.h \
      getColor.h \
      object.h \
      plane.h \
      shape.h \
      sphere.h \
      triple.h \
      window.h \
      world.h \
      polygon.h \
      light.h \
      rectangle.h \
      triangle.h

ray1: ${OBJS}
	@echo
	@echo Linking $@
	gcc -o $@ -Wall -g ${OBJS} -lm

${OBJS}: ${INCS} Makefile

test1:
	./ray1 world.txt output.ppm

test2:
	./ray1 world1.txt output.ppm

clean:
	rm -f *.o *.err ray1
