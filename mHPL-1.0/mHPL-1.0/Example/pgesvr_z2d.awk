
BEGIN {      nlines = 0;               


# Global Dictionary
   

	name_dictionary["tpzgesvr"] = "tpdgesvr";
	name_dictionary["complex"]  = "real";

	name_dictionary["pcgetrf"] = "psgetrf";

	name_dictionary["pzgemm"]  = "pdgemm";
	name_dictionary["pzgetrf"] = "pdgetrf";
	name_dictionary["pzgetrs"] = "pdgetrs";
	name_dictionary["pzelset"] = "pdelset";
	name_dictionary["pzgesvr"] = "pdgesvr";

	name_dictionary["pconvertz2c"] = "pconvertd2s";
	name_dictionary["pconvertc2z"] = "pconverts2d";

	name_dictionary["hpl_pcgesv"] = "hpl_psgesv";
	name_dictionary["hpl_pzgesv"] = "hpl_pdgesv";

	name_dictionary["hpl_cblacsinit"] = "hpl_sblacsinit";
       


}



{
	for (i in name_dictionary) {
		
		value = name_dictionary[i];

		gsub( i, value );
	};


	# ----------------------------------
        # save the input for post-processing
	# ----------------------------------

	nlines = nlines + 1;
	lines[ nlines ]  = $0;
}

END {
	for(i=1; i <= nlines; i += 1) {
		print lines[i];
	};
}
