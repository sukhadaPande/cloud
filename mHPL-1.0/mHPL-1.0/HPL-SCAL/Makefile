SHELL            = /bin/sh
#
arch             = UNKNOWN
#
include options
#
## Targets #############################################################
#
all              :
	-$(CP) options ./pshpldir
	-$(CP) options ./pdhpldir
	-$(CP) options ./pchpldir
	-$(CP) options ./pzhpldir
#
	-$(CP) LU.dat ./pshpldir
	-$(CP) LU.dat ./pdhpldir
	-$(CP) LU.dat ./pchpldir
	-$(CP) LU.dat ./pzhpldir
# 
	( $(CD) pshpldir;         $(MAKE) )
	( $(CD) pdhpldir;         $(MAKE) )
	( $(CD) pchpldir;         $(MAKE) )
	( $(CD) pzhpldir;         $(MAKE) )
#
clean            :
	( $(CD) pshpldir;         $(MAKE) clean)
	( $(CD) pdhpldir;         $(MAKE) clean)
	( $(CD) pchpldir;         $(MAKE) clean)
	( $(CD) pzhpldir;         $(MAKE) clean)
#
# ######################################################################
