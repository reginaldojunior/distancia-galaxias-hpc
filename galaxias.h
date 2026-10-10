/* galaxias.h - definicoes comuns da aula pratica (Pthreads)
 *
 * Problema: dado um conjunto de N galaxias em uma caixa cubica, calcular
 * estatisticas das distancias entre TODOS os pares: soma, media, minima,
 * maxima e um histograma. Custo: N*(N-1)/2 pares -> O(N^2).
 */
#ifndef GALAXIAS_H
#define GALAXIAS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

/* Windows (MSVC/MinGW) nao tem unistd.h, sysconf nem clock_gettime:
 * usa-se QueryPerformanceCounter e GetSystemInfo. */
#if defined(_WIN32)
#include <windows.h>
#else
#include <unistd.h>
#endif

#define LADO_CAIXA      1000.0                        /* lado da caixa (Mpc)          */
#define NUM_FAIXAS      10                            /* faixas do histograma         */
#define DIAGONAL_CAIXA  (LADO_CAIXA * 1.7320508075688772) /* diagonal da caixa         */

/* @ini:tipos */
typedef struct { double x, y, z; } Galaxia;            /* posicao no espaco (Mpc)      */

typedef struct {
    double   soma_distancias;     /* soma das distancias        */
    double   distancia_minima;    /* menor distancia             */
    double   distancia_maxima;    /* maior distancia             */
    uint64_t total_pares;         /* quantidade de pares         */
    uint64_t histograma[NUM_FAIXAS]; /* contagem por faixa      */
} Resultado;
/* @fim:tipos */

static inline double agora(void)
{
#if defined(_WIN32)
    LARGE_INTEGER frequencia, contador;
    QueryPerformanceFrequency(&frequencia);
    QueryPerformanceCounter(&contador);
    return (double)contador.QuadPart / (double)frequencia.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
}

/* Quantos nucleos a maquina tem (para o numero padrao de threads). */
static inline int num_processadores(void)
{
#if defined(_WIN32)
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    return info.dwNumberOfProcessors > 0 ? (int)info.dwNumberOfProcessors : 1;
#else
    long quantidade = sysconf(_SC_NPROCESSORS_ONLN);
    return quantidade > 0 ? (int)quantidade : 1;
#endif
}

/* Progresso no STDERR, para ficar claro que o programa esta calculando e nao
 * travado: sao N*(N-1)/2 pares e, num build sem otimizacao, isso leva dezenas
 * de segundos com a tela totalmente muda. Vai para o stderr de proposito: a
 * saida padrao continua limpa para o bench.sh e o verifica.sh lerem.
 * Nao chamar de dentro de uma thread (o estado e estatico). */
static inline void progresso(const char *rotulo, long parte, long total)
{
    static long proximo = 0;
    if (total <= 0 || parte < proximo) return;
    int ultimo = (parte + 1 >= total);          /* garante o 100% e a quebra de linha */
    fprintf(stderr, "\r%s: %3.0f%%%s", rotulo, 100.0 * (double)parte / (double)total,
            ultimo ? "\n" : "");
    fflush(stderr);
    proximo = parte + total / 50 + 1;           /* a cada ~2% */
    if (proximo > total) proximo = total;       /* nunca ultrapassa: senao o
                                                   ultimo passo seria pulado */
}

/* Gera posicoes pseudoaleatorias reprodutiveis (xorshift64*). */
static inline void gerar_galaxias(Galaxia *galaxias, long num_galaxias, uint64_t semente)
{
    uint64_t estado = semente;
    for (long indice = 0; indice < num_galaxias; indice++) {
        double coordenada[3];
        for (int eixo = 0; eixo < 3; eixo++) {
            estado ^= estado >> 12; estado ^= estado << 25; estado ^= estado >> 27;
            uint64_t sorteado = estado * 2685821657736338717ULL;
            coordenada[eixo] = LADO_CAIXA * (double)(sorteado >> 11) / 9007199254740992.0;
        }
        galaxias[indice].x = coordenada[0];
        galaxias[indice].y = coordenada[1];
        galaxias[indice].z = coordenada[2];
    }
}

/* @ini:distancia */
static inline double distancia(const Galaxia *galaxia_a, const Galaxia *galaxia_b)
{
    double delta_x = galaxia_a->x - galaxia_b->x;
    double delta_y = galaxia_a->y - galaxia_b->y;
    double delta_z = galaxia_a->z - galaxia_b->z;
    return sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
}
/* @fim:distancia */

static inline void resultado_inicializa(Resultado *resultado)
{
    resultado->soma_distancias = 0.0;
    resultado->distancia_minima = DIAGONAL_CAIXA;
    resultado->distancia_maxima = 0.0;
    resultado->total_pares = 0;
    for (int faixa = 0; faixa < NUM_FAIXAS; faixa++) resultado->histograma[faixa] = 0;
}

/* @ini:acumula */
static inline void resultado_acumula(Resultado *resultado, double distancia_par)
{
    int faixa = (int)(distancia_par / DIAGONAL_CAIXA * NUM_FAIXAS);
    if (faixa >= NUM_FAIXAS) faixa = NUM_FAIXAS - 1;
    resultado->soma_distancias += distancia_par;
    if (distancia_par < resultado->distancia_minima) resultado->distancia_minima = distancia_par;
    if (distancia_par > resultado->distancia_maxima) resultado->distancia_maxima = distancia_par;
    resultado->histograma[faixa]++;
    resultado->total_pares++;
}
/* @fim:acumula */

/* Combina dois resultados parciais (operacao de reducao). */
/* @ini:junta */
static inline void resultado_junta(Resultado *destino, const Resultado *parcial)
{
    destino->soma_distancias += parcial->soma_distancias;
    if (parcial->distancia_minima < destino->distancia_minima)
        destino->distancia_minima = parcial->distancia_minima;
    if (parcial->distancia_maxima > destino->distancia_maxima)
        destino->distancia_maxima = parcial->distancia_maxima;
    destino->total_pares += parcial->total_pares;
    for (int faixa = 0; faixa < NUM_FAIXAS; faixa++)
        destino->histograma[faixa] += parcial->histograma[faixa];
}
/* @fim:junta */

static inline void resultado_imprime(const Resultado *resultado, double tempo)
{
    printf("pares   = %llu\n", (unsigned long long)resultado->total_pares);
    printf("soma    = %.17g\n", resultado->soma_distancias);
    printf("media   = %.4f Mpc\n", resultado->total_pares ? resultado->soma_distancias / resultado->total_pares : 0.0);
    printf("dmin    = %.6f Mpc\n", resultado->distancia_minima);
    printf("dmax    = %.6f Mpc\n", resultado->distancia_maxima);
    printf("hist    =");
    for (int faixa = 0; faixa < NUM_FAIXAS; faixa++)
        printf(" %llu", (unsigned long long)resultado->histograma[faixa]);
    printf("\ntempo   = %.3f s\n", tempo);
    fflush(stdout);            /* aparecer imediatamente, mesmo se o processo
                                * for interrompido ou a saida for um pipe */
}

/* Confere se o numero de pares contados e N*(N-1)/2. */
static inline void resultado_verifica_pares(const Resultado *resultado, long num_galaxias)
{
    uint64_t pares_esperados = (uint64_t)num_galaxias * (uint64_t)(num_galaxias - 1) / 2;
    printf("esperado= %llu pares -> %s\n", (unsigned long long)pares_esperados,
           resultado->total_pares == pares_esperados ? "OK" : "ERRADO (corrida de dados!)");
    fflush(stdout);
}

/* Uso: ./programa [N] [NUM_THREADS] */
static inline void le_args(int argc, char **argv, long *num_galaxias, int *num_threads)
{
    *num_galaxias = (argc > 1) ? atol(argv[1]) : 20000;
    *num_threads   = (argc > 2) ? atoi(argv[2]) : num_processadores();
    if (*num_galaxias < 2)  *num_galaxias = 2;
    if (*num_threads < 1) *num_threads = 1;
}

#endif /* GALAXIAS_H */