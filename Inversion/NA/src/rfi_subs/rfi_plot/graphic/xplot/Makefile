#
#	Make discription file for xplot plot library 
#

include ../../../../makena.macros   # Edit this file to change location or
                                    # type of compilers (e.g. to gcc, g77)


xplot_subs: ps_comm.o ps_colr.o xpakw.o xpakb.o xtra.o cfopen.o libs

.f.o:
	$(F77) $*.f -c
.c.o:
	$(CC) $*.c -c

libs: 
	$(AR) rv ./lib/libpspak.a ps_comm.o ps_colr.o 
	$(AR) rv ./lib/libxpakw.a xpakw.o xtra.o cfopen.o
	$(AR) rv ./lib/libxpak.a xpakb.o xtra.o cfopen.o

