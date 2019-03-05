include options

all: pchpl pcscal

pchpl:
	$(F77) -DUSE_HPL=1 -fast -o pchpl \
	pchpl_driver.F \
         pcmatgen.f \
         pmatgeninc.f \
         pcluinfo.f \
         pcgetrrv.f \
         pclaschk.f \
         pclafchk.f \
	 $(MHPL_LIB)

pcscal:
	$(F77) -fast -o pcscal \
	pcludriver.f \
         pcmatgen.f \
         pmatgeninc.f \
         pcluinfo.f \
         pcgetrrv.f \
         pclaschk.f \
         pclafchk.f 

clean:
	-$(RM) *.o pchpl pcscal
