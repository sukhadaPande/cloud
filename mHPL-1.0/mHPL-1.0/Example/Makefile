include options

all: tpzgesvr tpdgesvr tpzgesvr_hpl tpdgesvr_hpl copy clean

copy:
	$(CP) tpzgesvr tpdgesvr tpzgesvr_hpl tpdgesvr_hpl \
	tpgesvr.dat $(EXAMPLE_BIN)

tpdgesvr.F: tpzgesvr.F pgesvr_z2d.awk
	$(AWK) -f pgesvr_z2d.awk < tpzgesvr.F > tpdgesvr.F

tpzgesvr_hpl: scalapack_mod.F pgesvr_mod.F \
	assert.f  profinit.f crandom_number.f tpzgesvr.F
	$(FC) $(FFLAGS) -DUSE_HPL=1 -o tpzgesvr_hpl \
	scalapack_mod.F pgesvr_mod.F \
		assert.f  profinit.f crandom_number.f \
		tpzgesvr.F \
		$(HPL_LIB) $(LIBS)

tpdgesvr_hpl: scalapack_mod.F pgesvr_mod.F \
	assert.f  profinit.f crandom_number.f tpdgesvr.F
	$(FC) $(FFLAGS) -DUSE_HPL=1 -o tpdgesvr_hpl \
	scalapack_mod.F pgesvr_mod.F \
		assert.f  profinit.f crandom_number.f \
		tpdgesvr.F \
		$(HPL_LIB) $(LIBS)

tpzgesvr: scalapack_mod.F pgesvr_mod.F \
	assert.f  profinit.f crandom_number.f tpzgesvr.F
	$(FC) $(FFLAGS) -o tpzgesvr \
	scalapack_mod.F pgesvr_mod.F \
		assert.f  profinit.f crandom_number.f \
		tpzgesvr.F \
		$(LIBS)

tpdgesvr: scalapack_mod.F pgesvr_mod.F \
	assert.f  profinit.f crandom_number.f tpdgesvr.F
	$(FC) $(FFLAGS) -o tpdgesvr \
	scalapack_mod.F pgesvr_mod.F \
		assert.f  profinit.f crandom_number.f \
		tpdgesvr.F \
		$(LIBS)
clean:
	$(RM) *.mod *.o tpzgesvr tpdgesvr tpzgesvr_hpl tpdgesvr_hpl

