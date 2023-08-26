MPICXX  := mpic++
CUDACXX := nvcc
IDIR    := include
LDIR    := lib
ODIR	:= build
SRCDIR := src
DOCDIR := documents

TARGET := DendriticGrowth_MPI.cpp
Nproc  := 1

CXXFLAGS := -I${IDIR}
LIBPATH  := -L/usr/local/cuda/lib64

_CUDADEPENDENCIES := General_CUDA_2D.cuh
CUDADEPENDENCIES  := ${patsubst %,${IDIR}/%,${_CUDADEPENDENCIES}}

_DEPENDENCIES := General.h General_MPI_2D.h
DEPENDENCIES  := ${patsubst %,${IDIR}/%,${_DEPENDENCIES}}

_OBJECTFILES := ${patsubst %.h,%.o,${_DEPENDENCIES}} ${patsubst %.cuh,%.o,${_CUDADEPENDENCIES}} 
_OBJECTFILES += ${patsubst %MPI.cpp,%CUDA.o,${TARGET}}
OBJECTFILES  := ${patsubst %,${ODIR}/%,${_OBJECTFILES}} ${patsubst %.cpp,${ODIR}/%.o,${TARGET}}

all: ${OBJECTFILES}
	${MPICXX} -std=c++17 -O3 $^ ${LIBPATH} -lcudart -o run

${ODIR}/%.o: %.cpp ${DEPENDENCIES}
	${MPICXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS}

${ODIR}/%.o: ${SRCDIR}/%.cpp ${DEPENDENCIES}
	${MPICXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS}

${ODIR}/General_CUDA_2D.o: ${SRCDIR}/General_CUDA_2D.cu ${CUDADEPENDENCIES}
	nvcc -c $< -o $@ ${CXXFLAGS}

${ODIR}/${patsubst %MPI.cpp,%CUDA.o,${TARGET}}: ${SRCDIR}/${patsubst %MPI.cpp,%CUDA.cu,${TARGET}} ${CUDADEPENDENCIES}
	nvcc -c $< -o $@ ${CXXFLAGS}

run:
	mpirun -np ${Nproc} ./run 2>/dev/null

.PHONY: clean run document

document:
	@xelatex -synctex=1 -interaction=nonstopmode ${DOCDIR}/Document.tex >> /dev/null
	@rm -rf ${DOCDIR}/*.aux ${DOCDIR}/*.log ${DOCDIR}/*.gz ${DOCDIR}/*.out
	@rm -rf *.aux *.log *.gz *.out
	@mv *.pdf ${DOCDIR}/

clean:
	rm -rf ${ODIR}/*.o ${SRCDIR}/*~ *~ run ${IDIR}/*~ Output/*.tec Output/*.vtk Output/*.visit
