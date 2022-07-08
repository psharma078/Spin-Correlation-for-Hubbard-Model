include /home/psharma/ITensor_fix/this_dir.mk
include /home/psharma/ITensor_fix/options.mk

#Define Flags ----------

TENSOR_HEADERS=$(PREFIX)/itensor/all.h
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L'$(ITENSOR_LIBDIR)' $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L'$(ITENSOR_LIBDIR)' $(ITENSOR_LIBGFLAGS)

#Rules ------------------

%.o: %.cc $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: Hubbard_model_rand_init_wf

debug: Hubbard_model_rand_init_wf

all: Hubbard_model_rand_init_wf

Hubbard_model_rand_init_wf: Hubbard_model_rand_init_wf.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) Hubbard_model_rand_init_wf.o -o Hubbard_model_rand_init_wf $(LIBFLAGS) -Wl,-rpath,$(ITENSOR_LIBDIR)

Hubbard_model_rand_init_wf-g: mkdebugdir .debug_objs/Hubbard_model_rand_init_wf.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/Hubbard_model_rand_init_wf.o -o Hubbard_model_rand_init_wf-g $(LIBGFLAGS) -Wl,-rpath,$(ITENSOR_LIBDIR)

mixedspin: mixedspin.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) mixedspin.o -o mixedspin $(LIBFLAGS)

mixedspin-g: mkdebugdir .debug_objs/mixedspin.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/mixedspin.o -o mixedspin-g $(LIBGFLAGS)

Hubbard_model_rand_init_wf_table: Hubbard_model_rand_init_wf_table.o $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) Hubbard_model_rand_init_wf_table.o -o Hubbard_model_rand_init_wf_table $(LIBFLAGS)

Hubbard_model_rand_init_wf_table-g: mkdebugdir .debug_objs/Hubbard_model_rand_init_wf_table.o $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCGFLAGS) .debug_objs/Hubbard_model_rand_init_wf_table.o -o Hubbard_model_rand_init_wf_table-g $(LIBGFLAGS)

mkdebugdir:
	mkdir -p .debug_objs

clean:
	@rm -fr *.o 
