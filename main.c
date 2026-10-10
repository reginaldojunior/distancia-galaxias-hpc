/* Versao sequencial: ponto de partida da aula. */
#include "galaxias.h"

/* @ini:kernel */
static Resultado calcula(const Galaxia *galaxias, long num_galaxias)
{
    Resultado resultado;
    resultado_inicializa(&resultado);
    for (long linha = 0; linha < num_galaxias - 1; linha++) {
        for (long coluna = linha + 1; coluna < num_galaxias; coluna++)
            resultado_acumula(&resultado, distancia(&galaxias[linha], &galaxias[coluna]));
        progresso("calculando", linha + 1, num_galaxias - 1);   /* so no stderr */
    }
    return resultado;
}
/* @fim:kernel */

/* @ini:main */
int main(int argc, char **argv)
{
    long num_galaxias; int num_threads;
    le_args(argc, argv, &num_galaxias, &num_threads);   /* num_threads e ignorado aqui */

    Galaxia *galaxias = malloc(num_galaxias * sizeof *galaxias);
    gerar_galaxias(galaxias, num_galaxias, 42);          /* fase 1: dados      */

    double tempo_inicio = agora();
    Resultado resultado = calcula(galaxias, num_galaxias);  /* fase 2: calculo */
    double tempo_fim = agora();

    printf("== sequencial: N = %ld ==\n", num_galaxias);
    resultado_imprime(&resultado, tempo_fim - tempo_inicio);
    resultado_verifica_pares(&resultado, num_galaxias);
    free(galaxias);
    return 0;
}
/* @fim:main */