# GRR20163049 Bruno Henrique Labres
# GRR20182981 Giovani Gurkevicz Marciniak

LIKWID = /home/soft/likwid
LIKWID_FLAGS = -I$(LIKWID)/include
LIKWID_LIBS = -L$(LIKWID)/lib

CFLAGS = -Wall -g $(LIKWID_FLAGS) -DLIKWID_PERFMON 
LDLIBS = -lm $(LIKWID_LIBS) -llikwid

objs = edp_lib.o pdeSolver.o

# regra default (primeira regra)
all: pdeSolver doc

# regras de ligacao
pdeSolver: $(objs)
	gcc -o pdeSolver $(objs) $(LDLIBS)
edp_lib.o: edp_lib.c edp_lib.h
	gcc -c edp_lib.c edp_lib.h $(CFLAGS) $(LDLIBS)

doc: dconfig
	doxygen $<

dconfig:
	doxygen -g dconfig

# remove arquivos temporários
clean:
	-rm $(objs)

# remove tudo o que não for o código-fonte
purge: clean
	-rm -r $(exec) doc
