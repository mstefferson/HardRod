
.PHONY : default
default :
	cd src; $(MAKE)
	
.PHONY : clean
clean :
	rm -f src/*.o 

force-build:

