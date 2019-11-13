# GRR20163049 Bruno Henrique Labres
# GRR20182981 Giovani Gurkevicz Marciniak

#CFLAGS = -Wall -g -std=c99  -O3 -mavx -march=native
CFLAGS = -Wall -g -std=c99 -O3 -mavx  -march=native -DLIKWID_PERFMON -I/home/soft/likwid/bin:/home/soft/likwid/sbin:${PATH}/include -L/home/soft/likwid/bin:/home/soft/likwid/sbin:${PATH}/lib -llikwid
LDLIBS = -lm

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
