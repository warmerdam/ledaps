#############################################################################
# !make
#
# makefile name: makefile (for linux)
#
##!END
#############################################################################
EXTRA  = -m32 -D_BSD_SOURCE

TARGET = lndpm

INC = -I. -I$(TIFFINC) -I$(GEOTIFF_INC)

LIB = -L$(GEOTIFF_LIB) -lgeotiff -L$(TIFFLIB) -ltiff -lm

OBJ =  lndpm.c

all: $(TARGET)

# Make the program
$(TARGET) : $(OBJ)
	$(CC) $(EXTRA) $(INC) $(OBJ) $(LIB)  -o $@

clean:
	rm -f *.o $(TARGET)

install:
	cp $(TARGET) $(BIN)

#******************* End of make file *******************************


