/* Pthreads v2 - CORRETA, mas lenta: um mutex protege o Resultado global. */
#include <pthread.h>
#include "galaxias.h"

typedef struct {
    const Galaxia *galaxias;
    long num_galaxias;
    int  id_thread;
    int  num_threads;
} Argumentos;

static Resultado resultado_global;
static pthread_mutex_t mutex_resultado = PTHREAD_MUTEX_INITIALIZER;

/* @ini:worker */
static void *trabalha(void *dados)
{
    const Argumentos *arg = dados;
    long inicio_linha = arg->id_thread * arg->num_galaxias / arg->num_threads;
    long fim_linha    = (arg->id_thread + 1) * arg->num_galaxias / arg->num_threads;

    for (long linha = inicio_linha; linha < fim_linha; linha++)
        for (long coluna = linha + 1; coluna < arg->num_galaxias; coluna++) {
            double distancia_par = distancia(&arg->galaxias[linha], &arg->galaxias[coluna]);
            pthread_mutex_lock(&mutex_resultado);   /* secao critica */
            resultado_acumula(&resultado_global, distancia_par);
            pthread_mutex_unlock(&mutex_resultado);
        }
    return NULL;
}
/* @fim:worker */

int main(int argc, char **argv)
{
    long num_galaxias; int num_threads;
    le_args(argc, argv, &num_galaxias, &num_threads);
    Galaxia *galaxias = malloc(num_galaxias * sizeof *galaxias);
    gerar_galaxias(galaxias, num_galaxias, 42);

    pthread_t *threads = malloc(num_threads * sizeof *threads);
    Argumentos *argumentos = malloc(num_threads * sizeof *argumentos);
    resultado_inicializa(&resultado_global);

    double tempo_inicio = agora();
    for (int id_thread = 0; id_thread < num_threads; id_thread++) {
        argumentos[id_thread] = (Argumentos){ galaxias, num_galaxias, id_thread, num_threads };
        pthread_create(&threads[id_thread], NULL, trabalha, &argumentos[id_thread]);
    }
    for (int id_thread = 0; id_thread < num_threads; id_thread++)
        pthread_join(threads[id_thread], NULL);
    double tempo_fim = agora();

    printf("== pthreads v2 (mutex): N = %ld, T = %d ==\n", num_galaxias, num_threads);
    resultado_imprime(&resultado_global, tempo_fim - tempo_inicio);
    resultado_verifica_pares(&resultado_global, num_galaxias);
    return 0;
}
