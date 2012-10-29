include ../Makefile.opt

NCFLAGS  = $(CFLAGS) $(LAC_FLAGS) $(INCDIR)
LIBS 	=	$(HDFLIB) $(TIFFLIB) $(BASELIB) 

TARGET1	= lndcal
OBJ1    = lndcal.o param.o space.o lut.o output.o cal.o util.o\
          names.o date.o myhdf.o mystring.o error.o input.o
INC1    = lndcal.h keyvalue.h param.h input.h lut.h output.h cal.h \
          names.h date.h myhdf.h mystring.h bool.h const.h error.h \
          space.h myproj.h myproj_const.h util.h

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
