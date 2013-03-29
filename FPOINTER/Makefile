

wrkdir/%.cu:		$(SRCDIR)/%.cu $(SRCDIR)/%.hh 
			@echo "#line 1 \"$<\"" > wrkdir/%.tmp 
			@cat wrkdir/%.tmp $< > $@
			@rm wrkdir/%.tmp
			@echo "Created $@"

wrkdir/CUDAglob.cu:	$(WRKFUNCTORLIST) $(HEADERLIST) 
			@rm -f $@
			@cat $(WRKFUNCTORLIST) > $@ 
