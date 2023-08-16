# Generated automatically from Makefile.in by configure.
#
# This is a Gromacs 3.0 template makefile for your own utility programs.
#
# Copy this file to whatever directory you are using for your own
# software and add more targets like the template one below.
#
# If you are using gmake it is relatively straightforward to add
# an include based on environment variables (like previous Gromacs versions)
# to select compiler flags and stuff automatically, but below it is static:
#

# Variables set by the configuration script:
#CC = icc
#LD = $(CC)

#XDRLIB = /usr/local/lib/libxdrfile.a
#CFLAGS = -O3 -w -funroll-all-loops -I/usr/local/include/xdrfile/ -lm
XDRLIB = /home/jhpeng/tools/xdrfile-1.1.4/lib/libxdrfile.a
#CFLAGS = -O3 -w -funroll-all-loops -I/usr/local/include/xdrfile/ -lm
CFLAGS = -O3 -w -funroll-all-loops -I/home/jhpeng/tools/xdrfile-1.1.4/include/ -lm
CC = gcc

# The real make targets - note that most make programs support
# the shortcut $^ instead of listing all object files a second
# time, but we cannot count on it...

unwrap:	unwrap.o smalloc.o
	$(CC) ${CFLAGS} -o $@ *.o $(XDRLIB)
clean:
		rm -f *.o
