SYSTEM  = IRIS4
#EXTRA	= -D$(SYSTEM) -lfastm
#EXTRA   = -D_BSD_SOURCE -Wall
EXTRA   = -m32 -D_BSD_SOURCE

INCDIR  = -I. -I$(TIFFINC) -I$(GEOTIFF_INC) -I$(HDFINC) -I$(HDFEOS_INC)
###INCDIR  = -I. -I$(HDFINC) -I$(HDFEOS_INC) -I$(TIFFINC) -I$(GEOTIFF_INC)
#NCFLAGS  = -mips3 -xansi -fullwarn -r10000 -g $(EXTRA) $(INCDIR)
#NCFLAGS  = $(CFLAGS) -r10000 -g $(EXTRA) $(INCDIR)
#NCFLAGS  = $(CFLAGS) -r10000 -O3 $(EXTRA) $(INCDIR)
#NCFLAGS  = $(CFLAGS) -g $(EXTRA) $(INCDIR)
#NCFLAGS  = $(CFLAGS) -O3 $(EXTRA) $(INCDIR)
#NCFLAGS  = $(CFLAGS) -O3  -DIRIS4 $(INCDIR)
NCFLAGS  = $(CFLAGS) -g $(INCDIR)

#EXLIB   = -lmalloc_ss
EXLIB	= -L$(GEOTIFF_LIB) -lgeotiff \
          -L$(TIFFLIB) -ltiff \
          -L$(HDFEOS_LIB) -lhdfeos -lgctp \
          -L$(HDFLIB) -lmfhdf -ldf -ljpeg -lz
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

TARGET1 = lndcsm
OBJ1    = lndcsm.o param.o input.o lut.o output.o csm.o space.o \
          names.o  myhdf.o mystring.o error.o tiff.o virbuf.o date.o util.o
INC1    = lndcsm.h keyvalue.h param.h input.h lut.h output.h csm.h names.h \
          date.h myhdf.h mystring.h bool.h const.h error.h tiff.h virbuf.h \
          util.h space.h myproj.h myproj_const.h


all: $(TARGET1)

x: $(TARGET1)

$(OBJ1): $(INC1)

$(TARGET1): $(OBJ1)
	$(CC) $(EXTRA) -o $(TARGET1) $(OBJ1) $(LOADLIB)
#	/usr/pure/purify -m71-engine $(CC) -o $(TARGET1) $(OBJ1) $(LOADLIB) 

clean:
	rm -f *.o $(TARGET1)

install:
	cp $(TARGET1) $(BIN)

#
# Rules
#
.c.o:
	$(CC) $(EXTRA) $(NCFLAGS) -c $< -o $@
