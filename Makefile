GXX = g++
SLD =
DLD = -lpthread
OBJ = qtar_opts.o qtar_subs.o qtar_core.o qtar.o
VPATH = src

qtar : ${OBJ}
	${GXX} -o qtar src/sw/ssw_cpp.cpp src/sw/ssw.c ${OBJ} ${DLD} ${SLD}
.PHONY : clean
clean :
	rm -f ${OBJ}

qtar_opts.o:qtar_opts.h
qtar_subs.o:qtar_subs.h
qtar_core.o:qtar_core.h
qtar.o:qtar.h
