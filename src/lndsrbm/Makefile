FC      = gfortran
EXTRA   = -m32 -D_BSD_SOURCE
INCDIR  = -I. -I$(HDFINC)
EXLIB   = -L$(HDFLIB) -lmfhdf -ldf -ljpeg -lz -lm
EXE     = cmrbv1.0 compadjn comptemp xy2geo geo2xy SDSreader3.0
GEOLOC_DEPEND = LS_geoloc.o LS_geoloc_driver.c

all : $(EXE)

cmrbv1.0 : cmrbv1.0.f
	$(FC) $(EXTRA) -g -o $@ $? $(INCDIR) $(EXLIB)

compadjn : compadjn.f
	$(FC) $(EXTRA) -g -o $@ $?

comptemp : comptemp.f
	$(FC) $(EXTRA) -g -o $@ $?

LS_geoloc.o : LS_geoloc.c
	$(CC) $(EXTRA) -O0 -o $@ -c $?

xy2geo : $(GEOLOC_DEPEND)
	$(CC) $(EXTRA) -DINV -O0 -o $@ $(GEOLOC_DEPEND) $(INCDIR) $(EXLIB)

geo2xy : $(GEOLOC_DEPEND)
	$(CC) $(EXTRA) -O0 -o $@ $(GEOLOC_DEPEND) $(INCDIR) $(EXLIB)

SDSreader3.0 : SDSreader3.0.c
	$(CC) $(EXTRA) -O0 -o $@ $? $(INCDIR) $(EXLIB)

c.o : 
	$(CC) $(EXTRA) -O0 -o $@ -c $?

clean :
	rm -f *.o
	rm -f $(EXE)

install:
	cp cmrbv1.0 $(BIN)
	cp compadjn $(BIN)
	cp comptemp $(BIN)
	cp xy2geo $(BIN) 
	cp geo2xy $(BIN)
	cp SDSreader3.0 $(BIN)
	cp lndsrbm.ksh $(BIN)

