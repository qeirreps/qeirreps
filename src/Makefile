# Copyright (c) Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe, 2020
FC     = ifort
LFLAGS = -mkl

.SUFFIXES: .f90

.f90.o:
	$(FC) -c $<

TARGET = qeirreps.x
OBJS = qeirreps.o

${TARGET} : $(OBJS)
	$(FC) $(LFLAGS) -o $@ $(OBJS) 

clean:
	${RM} ${TARGET} ${OBJS}
