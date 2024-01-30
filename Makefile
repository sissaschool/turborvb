$(info ===================================================)
$(info )
$(info TurboRVB Legacy build system)
$(info )
$(info ===================================================)
$(info )

# Test if make.inc exists

make_inc_c := e
make_txt_c := e

ifeq ($(wildcard make.inc),)
$(info )
$(warning make.inc not found. )

make_inc_c := $(shell read -p "Do you want to create make.inc? [y/n]: " ans; \
                  if [ "$$ans" = "y" ]; then \
                      cp devel_tools/configs/make.inc.example.gcc make.inc; \
                      echo "c"; \
                  else \
                      echo "n"; \
                  fi)

ifeq ($(make_inc_c) ,c)
    $(info )	
    $(info Please select an example: )
    $(info )
    $(info 1 GNU Fortran Compiler )
    $(info 2 GNU Fortran Compiler with MPI) 
    $(info )
    $(shell read -p "Select an example [1-2]: " ans; \
        if [ "$$ans" = "1" ]; then \
            cp devel_tools/make.inc.examples/make.inc.example.gcc make.inc; \
        elif [ "$$ans" = "2" ]; then \
            cp devel_tools/configs/make.inc.example.gccmpi make.inc; \
        fi)

endif

endif

ifeq ($(wildcard make.txt),)
$(info )
$(warning make.txt not found. )

make_txt_c := $(shell read -p "Do you want to create make.txt? [y/n]: " ans; \
                  if [ "$$ans" = "y" ]; then \
                      cp devel_tools/configs/make.txt.example make.txt; \
                      echo "c"; \
                  else \
                      echo "n"; \
                  fi)

endif

$(info )

ifeq ($(make_inc_c),c)
$(info make.inc created. Please edit it.)
endif

ifeq ($(make_txt_c),c)
$(info make.txt created. Please edit it.)
endif

ifeq ($(make_inc_c),n)
$(info make.inc not created. Exiting.)
endif

ifeq ($(make_inc_c),n)
$(info make.inc not created. Exiting.)
endif

ifneq ($(make_inc_c),e)
$(info )
$(error Exiting.)
endif

ifneq ($(make_txt_c),e)
$(info )
$(error Exiting.)
endif

include make.inc

$(info )
$(info Using make.inc found in $(CURDIR))
$(info Building in $(BUILD_DIR))
$(info )
$(info Fortran compiler: $(FC))
$(info Fortran compiler flags: $(FCFLAGS))
$(info Fortran passive flags: $(FCFLAGS_PASSIVE))
$(info Fortran aggressive flags: $(FCFLAGS_AGGRESSIVE))
$(info )
$(info C compiler: $(CC))
$(info C compiler flags: $(CFLAGS))
$(info )
$(info Linker: $(LD))
$(info Linker flags: $(FLINK))
$(info )
$(info Libraries: $(LINK_LIBS))
$(info )
$(info ===================================================)


all:
	make -C src/m_pfapack -f Makefile
	make -C src/m_common -f Makefile
	make -C src/d_qlapack -f Makefile
	make -C src/c_adjoint_forward -f Makefile
	make -C src/c_adjoint_backward -f Makefile
	make -C src/b_complex -f Makefile
	make -C src/a_turborvb -f Makefile
	
	@echo " "
	@echo " DONE - Enjoy TurboRVB! "
	@echo " "

clean:
	rm -rf $(BUILD_DIR)

