CD=cd
RM=rm
CP=cp
MKDIR=mkdir
SHELL=bash

all:
	$(MKDIR) -p lib bin include
	cd src && $(MAKE)

clean:
	cd src && $(MAKE) clean

realclean:
	cd src && $(MAKE) clean
	$(RM) -f bin/*
	$(RM) -f lib/*
	$(RM) -f include/*

