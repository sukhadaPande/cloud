include options

all: pshpl psscal

pshpl:
	$(F77) -DUSE_HPL=1 -fast -o pshpl \
	pshpl_driver.F \
         psmatgen.f \
         pmatgeninc.f \
         psluinfo.f \
         psgetrrv.f \
         pslaschk.f \
         pslafchk.f \
	 $(MHPL_LIB)

psscal:
	$(F77) -fast -o psscal \
	psludriver.f \
         psmatgen.f \
         pmatgeninc.f \
         psluinfo.f \
         psgetrrv.f \
         pslaschk.f \
         pslafchk.f 

clean:
	-$(RM) *.o pshpl psscal
