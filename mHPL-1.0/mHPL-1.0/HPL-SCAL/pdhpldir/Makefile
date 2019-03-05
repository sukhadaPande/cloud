include options

all: pdhpl pdscal

pdhpl:
	$(F77) -DUSE_HPL=1 -fast -o pdhpl \
	pdhpl_driver.F \
         pdmatgen.f \
         pmatgeninc.f \
         pdluinfo.f \
         pdgetrrv.f \
         pdlaschk.f \
         pdlafchk.f \
	 $(MHPL_LIB)

pdscal:
	$(F77) -fast -o pdscal \
	pdludriver.f \
         pdmatgen.f \
         pmatgeninc.f \
         pdluinfo.f \
         pdgetrrv.f \
         pdlaschk.f \
         pdlafchk.f 

clean:
	-$(RM) *.o pdhpl pdscal
