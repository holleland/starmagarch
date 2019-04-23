all:	STARMAGARCH.so

STARMAGARCH.so:	STARMAGARCH.cpp
    Rscript --vanilla -e "TMB::compile('STARMAGARCH.cpp','-O0 -g')"

clean:	rm -rf *o
