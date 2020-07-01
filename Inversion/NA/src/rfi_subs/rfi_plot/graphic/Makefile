#
#	Make discription file for Shibutani plot library 
#

#FFLAGS = -O4
#F77 = /opt/local/SUNWspro/bin/f77

include ../../../makena.macros   # Edit this file to change location or
                                 # type of compilers (e.g. to gcc, g77)


all: glib

glib: arc.o arcn.o axis.o axs.o ccirc1.o cellsize.o circ1.o circle.o \
      csymbl.o dashln.o fcell.o fcellg.o frame.o newpage.o newpen.o \
      number.o origin.o penw.o plot.o plote.o plots.o plus.o rot.o \
      rotate.o scale.o sqr.o star.o symbol.o tri.o xnumber.o \
      xplot.o xsymbol.o 

.f.o:
	$(F77) $*.f -c
	$(AR) rv lib/glib.a $*.o

