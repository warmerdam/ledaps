SYSTEM  = IRIS4
#EXTRA	= -D$(SYSTEM) -lfastm
#EXTRA   = -DLINUX -D_BSD_SOURCE -Wall
EXTRA   = -DLINUX -D_BSD_SOURCE -m32

INCDIR  = -I. -I$(HDFINC) -I$(HDFEOS_INC) -I$(TIFFINC) -I$(GEOTIFF_INC)
#NCFLAGS  = -mips3 -xansi -fullwarn -r10000 -g $(EXTRA) $(INCDIR)
#NCFLAGS  = $(CFLAGS) -r10000 -g $(EXTRA) $(INCDIR)
#NCFLAGS  = $(CFLAGS) -r10000 -O3 $(EXTRA) $(INCDIR)
NCFLAGS  = $(CFLAGS) -O $(EXTRA) $(INCDIR)
#NCFLAGS  = $(CFLAGS) $(EXTRA) $(INCDIR)

#EXLIB   = -lmalloc_ss
EXLIB	= -L$(HDFEOS_LIB) -lhdfeos -lgctp \
          -L$(HDFLIB) -lmfhdf -ldf -lz -ljpeg
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

TARGET1	= lndsr
#EXTTARGET = sixsV1.0B

OBJ1    = lndsr.o param.o input.o prwv_input.o lut.o output.o sr.o ar.o \
          names.o date.o myhdf.o mystring.o error.o space.o grib.o read_grib_tools.o \
          CHAND.o CSALBR.o sixs_runs.o clouds.o

INC1    = lndsr.h keyvalue.h param.h input.h prwv_input.h lut.h output.h sr.h ar.h \
          names.h date.h myhdf.h mystring.h bool.h const.h error.h space.h grib.h read_grib_tools.h \
          myproj.h myproj_const.h

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
#	cp $(EXTTARGET) $(BIN)

#
# Rules
#
.c.o:
	$(CC) $(EXTRA) $(NCFLAGS) -c $< -o $@
.f.o:
	gfortran -m32 -c $< -o $@
