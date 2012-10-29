include ../Makefile.opt

NCFLAGS  = $(CFLAGS) $(LAC_FLAGS) $(INCDIR)
LIBS 	=	$(HDFLIB) $(BASELIB) 

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
	$(CC) $(NCFLAGS) -o $(TARGET1) $(OBJ1) $(LIBS)
#	/usr/pure/purify -m71-engine $(CC) -o $(TARGET1) $(OBJ1) $(LOADLIB) 

clean:
	rm -f *.o $(TARGET1)

install:	$(TARGET1)
	cp $(TARGET1) $(BIN)

#
# Rules
#
.c.o:
	$(CC) $(NCFLAGS) -c $< -o $@
