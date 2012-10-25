#############################################################################
# !make
#
# makefile name: lndapp.mk (for linux)
#
##!END
#############################################################################
CC = gcc -m32
CFLAGS = -w

TARGET = lndapp

INC = -I. -I$(HDFINC) -I$(HDFEOS_INC)

LIB = -L$(HDFEOS_LIB) -lhdfeos -lgctp\
      -L$(HDFLIB) -lmfhdf -ldf -ljpeg -lz -lm

OBJ = lndapp.o

all: $(TARGET)

$(OBJ) : $(INC_FILES)

# Make the process
$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET) $(LIB)

clean:
	rm -f *.o $(TARGET)

install:
	cp $(TARGET) $(BIN)

#
# Rules
#
.c.o: $(INC_FILES)
	$(CC) $(CFLAGS) $(ADD_CFLAGS) $(INC) -c $< -o $@

#******************* End of make file *******************************


