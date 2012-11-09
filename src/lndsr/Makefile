include ../Makefile.opt

NCFLAGS  = $(CFLAGS) $(LAC_FLAGS) $(INCDIR) \
	-DSAVE_6S_RESULTS -DSAVE_AR_RESULTS
LIBS = $(HDFLIB) $(BASELIB)

TARGET1	= lndsr
#EXTTARGET = sixsV1.0B

CORE_OBJ    = param.o input.o prwv_input.o lut.o output.o sr.o ar.o \
          names.o date.o myhdf.o mystring.o error.o space.o grib.o read_grib_tools.o \
          CHAND.o CSALBR.o sixs_runs.o clouds.o lndsr_util.o

INC1    = lndsr.h keyvalue.h param.h input.h prwv_input.h lut.h output.h sr.h ar.h \
          names.h date.h myhdf.h mystring.h bool.h const.h error.h space.h grib.h read_grib_tools.h \
          myproj.h myproj_const.h

OBJ_EE  =	ee_lndsr_test.o ee_lndsr.o ee_gridcell.o ee_inputmeta.o $(CORE_OBJ)
OBJ1	=	lndsr.o $(CORE_OBJ)

all: $(TARGET1) ee_lndsr

x: $(TARGET1)

$(OBJ_EE) $(OBJ1): $(INC1) ee_lndsr.h Makefile ../Makefile.opt

$(TARGET1): $(OBJ1)
	$(CC) $(NCFLAGS) -o $(TARGET1) $(OBJ1) $(LIBS)
#	/usr/pure/purify -m71-engine $(CC) -o $(TARGET1) $(OBJ1) $(LIBS) 

ee_lndsr: $(OBJ_EE)
	$(CC) $(NCFLAGS) -o ee_lndsr $(OBJ_EE) $(LIBS)

clean:
	rm -f *.o $(TARGET1) ee_lndsr

install: $(TARGET1) ee_lndsr
	cp $(TARGET1) $(BIN)
#	cp $(EXTTARGET) $(BIN)
	cp ee_lndsr $(BIN)

#
# Rules
#
.c.o:
	$(CC) $(NCFLAGS) -c $< -o $@
.cc.o:
	$(CXX) $(NCFLAGS) -c $< -o $@
.f.o:
	gfortran $(LAF_FLAGS) -c $< -o $@
