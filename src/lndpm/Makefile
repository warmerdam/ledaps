#############################################################################
# !make
#
# makefile name: makefile (for linux)
#
##!END
#############################################################################

include ../Makefile.opt

TARGET = lndpm

NCFLAGS  = $(CFLAGS) $(LAC_FLAGS) $(INCDIR)
LIBS = $(HDFLIB) $(TIFFLIB) $(BASELIB)

OBJ =  lndpm.c

all: $(TARGET)

# Make the program
$(TARGET) : $(OBJ)
	$(CC) $(NCFLAGS) $(OBJ) $(LIBS)  -o $@

clean:
	rm -f *.o $(TARGET)

install:	$(TARGET)
	cp $(TARGET) $(BIN)

#******************* End of make file *******************************


