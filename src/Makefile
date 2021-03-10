FC = gfortran
FFLAG = -O2
AR = ar r
RM = rm
INST = ../bin
LPATH = lib
PROGS = qvib.x qfit.x
#
.PHONY:	all
all:	$(PROGS)
	mv *.x $(INST)
#
.PHONY:	clean
clean:	
	for f in $(PROGS) ; do rm $(INST)/$${f} ; done
	rm *.o lib/libc8v.a
	rm Makefile
#
qvib.x: qvib.f
	$(FC) $(FFLAG) -o qvib.x qvib.f $(LPATH)/qher.f $(LPATH)/qmat.f $(LPATH)/diag.f
qfit.x: qfit.f
	$(FC) $(FFLAG) -o qfit.x qfit.f
