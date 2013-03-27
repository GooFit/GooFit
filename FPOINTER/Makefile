

build/%.cu:		$(SRCDIR)/%.cu $(SRCDIR)/%.hh 
	@echo "#line 0 \"$<\"" > build/%.tmp 
	@cat build/%.tmp $< > $@
	@rm build/%.tmp
	@echo "Created $@"

build/CUDAglob.cu:	$(WRKFUNCTORLIST) $(HEADERLIST) 
	@rm -f $@
	@cat $(WRKFUNCTORLIST) > $@ 
