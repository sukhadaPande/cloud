include options

all: pzhpl pzscal

pzhpl:
	$(F77) -DUSE_HPL=1 -fast -o pzhpl \
	pzhpl_driver.F \
         pzmatgen.f \
         pmatgeninc.f \
         pzluinfo.f \
         pzgetrrv.f \
         pzlaschk.f \
         pzlafchk.f \
	 $(MHPL_LIB)

pzscal:
	$(F77)  -fast -o pzscal \
	pzludriver.f \
         pzmatgen.f \
         pmatgeninc.f \
         pzluinfo.f \
         pzgetrrv.f \
         pzlaschk.f \
         pzlafchk.f 

clean:
	-$(RM) *.o pzhpl pzscal
