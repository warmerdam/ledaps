include ../Makefile.opt

NCFLAGS  = $(CFLAGS) $(LAC_FLAGS) $(INCDIR)
LIBS = $(HDFLIB) $(BASELIB)

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

$(OBJ1): $(INC1) Makefile ../Makefile.opt

$(TARGET1): $(OBJ1)
	$(CC) $(NCFLAGS) -o $(TARGET1) $(OBJ1) $(LIBS)
#	/usr/pure/purify -m71-engine $(CC) -o $(TARGET1) $(OBJ1) $(LIBS) 

clean:
	rm -f *.o $(TARGET1)

install: $(TARGET1)
	cp $(TARGET1) $(BIN)
#	cp $(EXTTARGET) $(BIN)

#
# Rules
#
.c.o:
	$(CC) $(NCFLAGS) -c $< -o $@
.f.o:
	gfortran $(LAF_FLAGS) -c $< -o $@
