/* Pthreads v6 - MPI: memoria DISTRIBUIDA entre processos.
 *
 * Diferencas em relacao a v3/v4/v5 (pthreads):
 *   - nao existe memoria compartilhada: cada processo (rank) tem sua copia;
 *   - os dados NAO sao enviados: cada rank gera o mesmo conjunto (mesma semente),
 *     que em um caso real viria de um MPI_Scatter ou MPI_Bcast;
 *   - a barreira e um MPI_Barrier (equivalente ao pthread_join das threads);
 *   - a reducao final e um MPI_Gather seguido de resultado_junta no rank 0,
 *     porque Resultado e uma struct com tipos diferentes (soma/min/max/histograma)
 *     e nao pode ser reduzida direto pelo MPI_Reduce. A alternativa seria
 *     decompor em varios MPI_Reduce (MPI_SUM para soma e pares,
 *     MPI_MIN/MPI_MAX para as distancias, MPI_SUM para o histograma).
 *   - o tempo e medido com MPI_Wtime (o relogio proprio do MPI).
 *
 * Uso:  mpirun -np T ./galaxias_v6_mpi N
 * Compilar com particao em BLOCOS (igual a v3, desbalanceada):
 *        mpicc -O2 -Wall -Wextra -DPARTICAO_BLOCOS -o v6 galaxias_v6_mpi.c -lm
 */
#include <mpi.h>
#include "galaxias.h"

/* Cada rank calcula as estatisticas das linhas que lhe cabem. */
static void processa_linhas(const Galaxia *galaxias, long num_galaxias,
                            int rank, int tamanho, Resultado *parcial)
{
    resultado_inicializa(parcial);

#ifdef PARTICAO_BLOCOS
    long inicio_linha = rank * num_galaxias / tamanho;         /* bloco contiguo */
    long fim_linha    = (rank + 1) * num_galaxias / tamanho;
    for (long linha = inicio_linha; linha < fim_linha; linha++)
#else
    for (long linha = rank; linha < num_galaxias; linha += tamanho)   /* ciclico */
#endif
    {
        for (long coluna = linha + 1; coluna < num_galaxias; coluna++)
            resultado_acumula(parcial, distancia(&galaxias[linha], &galaxias[coluna]));
    }
}

/* @ini:main */
int main(int argc, char **argv)
{
    int rank, tamanho;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &tamanho);

    long num_galaxias; int num_threads_ignorado;   /* o paralelismo vem do -np */
    le_args(argc, argv, &num_galaxias, &num_threads_ignorado);

    Galaxia *galaxias = malloc(num_galaxias * sizeof *galaxias);
    gerar_galaxias(galaxias, num_galaxias, 42);    /* todos geram a MESMA caixa */

    MPI_Barrier(MPI_COMM_WORLD);                  /* ninguem mede enquanto outro cria dados */

    double tempo_inicio = MPI_Wtime();
    Resultado resultado_local;
    processa_linhas(galaxias, num_galaxias, rank, tamanho, &resultado_local);
    double tempo_rank = MPI_Wtime() - tempo_inicio;

    /* o rank 0 recebe os resultados parciais e os tempos de cada rank */
    Resultado *parciais = NULL;
    double *tempos_ranks = NULL;
    if (rank == 0) {
        parciais     = malloc(tamanho * sizeof *parciais);
        tempos_ranks = malloc(tamanho * sizeof *tempos_ranks);
    }
    MPI_Gather(&resultado_local, sizeof(Resultado), MPI_BYTE,
               parciais, sizeof(Resultado), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Gather(&tempo_rank, 1, MPI_DOUBLE,
               tempos_ranks, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);                  /* a barreira final: o rank 0 espera todos */
    double tempo_fim = MPI_Wtime();

    if (rank == 0) {
        Resultado resultado_total;                /* reducao identica a das pthreads */
        resultado_inicializa(&resultado_total);
        for (int r = 0; r < tamanho; r++)
            resultado_junta(&resultado_total, &parciais[r]);

        printf("== MPI v6 (%s): N = %ld, ranks = %d ==\n",
#ifdef PARTICAO_BLOCOS
               "blocos",
#else
               "ciclico",
#endif
               num_galaxias, tamanho);
        resultado_imprime(&resultado_total, tempo_fim - tempo_inicio);
        resultado_verifica_pares(&resultado_total, num_galaxias);
        for (int r = 0; r < tamanho; r++)
            printf("rank %d: pares = %10llu (%5.1f%%)  tempo = %.3f s\n", r,
                   (unsigned long long)parciais[r].total_pares,
                   100.0 * parciais[r].total_pares / resultado_total.total_pares,
                   tempos_ranks[r]);
        free(parciais); free(tempos_ranks);
    }

    free(galaxias);
    MPI_Finalize();
    return 0;
}
/* @fim:main */
