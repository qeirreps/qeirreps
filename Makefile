# Copyright (c) Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe
FC     = mpif90
FFLAGS = -qopenmp 
LFLAGS = -mkl -qopenmp 

.SUFFIXES: .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

TARGET = qeirreps.x
OBJS = qeirreps.o

${TARGET} : $(OBJS)
	$(FC) $(LFLAGS) -o $@ $(OBJS) ${LFLAGS}

clean:
	${RM} ${TARGET} ${OBJS}
