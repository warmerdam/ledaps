
MODLIST = lndpm lndcal lndcsm lndsr lndapp lndsrbm 6sV-1.0B

all:
	@for i in $(MODLIST); do \
	(cd $$i; make) \
	done

install:
	@for i in $(MODLIST); do \
	(cd $$i; make install) \
	done

clean:
	@for i in $(MODLIST); do \
	(cd $$i; make clean > /dev/null) \
	done
