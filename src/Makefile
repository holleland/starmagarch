all:	STARMAGARCH.so

STARMAGARCH.so:	STARMAGARCH.cpp
    Rscript --vanilla -e "TMB::compile('STARMAGARCH.cpp','-00 -g')"

clean:	rm -rf *o
