
 SUBDIRS = prog user

MAKESUBDIRS = @for dir in $(SUBDIRS); do \
			echo Making $@ in $$dir; \
			$(MAKE)	-C $$dir $@ || exit 1 ; done


all postscript:		
	$(MAKESUBDIRS)

pdf:
	$(MAKESUBDIRS)

html:
	$(MAKESUBDIRS)

clean:		
	$(MAKESUBDIRS)

realclean:	
	$(MAKESUBDIRS)
