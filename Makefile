include make.inc

all:
	make -C src/m_pfapack -f Makefile
	make -C src/m_common -f Makefile
	make -C src/d_qlapack -f Makefile
	make -C src/c_adjoint_forward -f Makefile
	make -C src/c_adjoint_backward -f Makefile
	make -C src/b_complex -f Makefile
