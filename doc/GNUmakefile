
 SUBDIRS = prog user nwchem_func

MAKESUBDIRS = @for dir in $(SUBDIRS); do \
			echo Making $@ in $$dir; \
			$(MAKE)	-C $$dir $@ || exit 1 ; done


all postscript:		
	$(MAKESUBDIRS)

html:
	$(MAKESUBDIRS)

clean:		
	$(MAKESUBDIRS)

realclean:	
	$(MAKESUBDIRS)
