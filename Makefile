CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -pthread
LDLIBS  = -lm
PROGS   = galaxias_seq galaxias_v1_corrida galaxias_v2_mutex \
          galaxias_v3_bloco galaxias_v4_ciclico galaxias_v5_dinamico

# O clang da Apple nao aceita -fopenmp: no macOS usa-se o libomp do Homebrew
# (prefixando o pragma com -Xpreprocessor).  Em Linux, -fopenmp direto.
LIBOMP_DIR ?= /opt/homebrew/opt/libomp
UNAME_S    := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
OMP_FLAGS = -Xpreprocessor -fopenmp -I$(LIBOMP_DIR)/include -L$(LIBOMP_DIR)/lib -lomp
else
OMP_FLAGS = -fopenmp
endif

MPICC     ?= mpicc
NVCC      ?= nvcc
CUDA_ARCH ?= sm_60          # atomicAdd(double*) exige CC 6.0+

all: $(PROGS) galaxias_v7_openmp galaxias_v6_mpi

$(PROGS): %: %.c galaxias.h pth_comum.h
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

# OpenMP: as mesmas threads da v3-v5, com o laco em pragma.
# Troque a particao recompilando, por exemplo para reproduzir a v4 (ciclica):
#   make galaxias_v7_openmp OMP_EXTRA="-DTIPO_ESCALA=omp_sched_static -DTAMANHO_BLOCO=1"
galaxias_v7_openmp: galaxias_v7_openmp.c galaxias.h
	$(CC) -O2 -Wall -Wextra $(OMP_FLAGS) $(OMP_EXTRA) -o $@ $< $(LDLIBS)

# MPI: memoria distribuida. Nao entra em PROGS porque nem toda maquina tem MPI:
#   make galaxias_v6_mpi
galaxias_v6_mpi: galaxias_v6_mpi.c galaxias.h
	$(MPICC) -O2 -Wall -Wextra -o $@ $< $(LDLIBS)

# CUDA: so compila em maquina com GPU NVIDIA (nao tem nvcc em Apple Silicon),
# por isso tambem fica fora de PROGS:
#   make galaxias_v8_cuda CUDA_ARCH=sm_75
galaxias_v8_cuda: galaxias_v8_cuda.cu galaxias.h
	$(NVCC) -O2 -arch=$(CUDA_ARCH) -o $@ $< $(LDLIBS)

# Versao para o gprof (sem inline, para os nomes das funcoes aparecerem)
perfil: galaxias_seq.c galaxias.h
	$(CC) -O2 -fno-inline -pg -o galaxias_seq_pg galaxias_seq.c $(LDLIBS)
	./galaxias_seq_pg 15000
	gprof -b -p ./galaxias_seq_pg gmon.out | head -12

clean:
	rm -f $(PROGS) galaxias_v7_openmp galaxias_v6_mpi galaxias_v8_cuda \
	      galaxias_seq_pg gmon.out
