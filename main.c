/*******************************************************************************
 *
 * DESCRIÇÃO GERAL:
 * Implementação de um sistema para encontrar o trajeto mais rápido entre
 * duas estações do metrô de Paris. O sistema utiliza o Algoritmo A* (busca informada),
 * calculando o tempo real de deslocamento nos trilhos (Tabela 2), somando o tempo
 * de baldeação (4 minutos por troca de linha) e utilizando como heurística admissível
 * a distância em linha reta (Tabela 1). Também permite comparar com o Algoritmo de
 * Dijkstra (busca cega com h = 0).
 *
 * ORGANIZAÇÃO DO CÓDIGO (DIVIDIDO EM 5 BLOCOS PARA APRESENTAÇÃO EM EQUIPE):
 * - BLOCO 1: Definição de Estruturas, Constantes e Dados do Grafo (Tabelas 1, 2 e Linhas)
 * - BLOCO 2: Cálculos de Tempo, Heurística Admissível e Transições com Baldeação
 * - BLOCO 3: Mecanismo do Algoritmo A* e Dijkstra (Seleção do Menor f e Expansão)
 * - BLOCO 4: Reconstrução do Percurso e Relatório Detalhado da Rota
 * - BLOCO 5: Interface com Usuário, Demonstração, Respostas Teóricas e Função Main
 *
 * INTEGRANTES DA EQUIPE:
 * 1. Bloco 1: Estruturas de Dados e Modelagem do Grafo
 * 2. Bloco 2: Custos de Deslocamento e Função Heurística
 * 3. Bloco 3: Algoritmo de Busca A* e Dijkstra
 * 4. Bloco 4: Reconstrução e Exibição do Percurso
 * 5. Bloco 5: Interface, Casos de Teste e Questões Teóricas
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*******************************************************************************
 * =============================================================================
 * BLOCO 1: ESTRUTURAS DE DADOS, CONSTANTES E MODELAGEM DO GRAFO
 * =============================================================================
 * Este bloco define as constantes fundamentais do sistema, os identificadores
 * das linhas de metrô, as tabelas com os dados geográficos e de malha viária
 * fornecidos na especificação do projeto, bem como a estrutura que armazena os
 * estados da busca.
 *******************************************************************************/

/* Quantidade total de estações na malha considerada (E1 a E14) */
#define NUM_ESTACOES 14

/* Quantidade total de linhas metroviárias */
#define NUM_LINHAS 4

/* Inclui o estado 0 (quando nenhuma linha foi utilizada na partida) + 4 linhas */
#define NUM_ESTADOS_LINHA (NUM_LINHAS + 1)

/* Velocidade média operacional dos trens em km/h */
#define VELOCIDADE_KM_H 30.0

/* Tempo de espera e caminhada para realizar baldeação (troca de linha) em minutos */
#define TEMPO_BALDEACAO_MIN 4.0

/* Representação de valor infinito para inicialização de custos */
#define VALOR_INFINITO 999999.0

/* Códigos identificadores das cores das linhas */
#define LINHA_NENHUMA  0
#define LINHA_AZUL     1
#define LINHA_AMARELA  2
#define LINHA_VERMELHA 3
#define LINHA_VERDE    4

/*
 * TABELA 1: DISTÂNCIAS DIRETAS (EM LINHA RETA) ENTRE ESTAÇÕES (em km)
 * Utilizada exclusivamente para o cálculo da função heurística h(n).
 * Matriz simétrica 14x14 contendo os valores exatos da Tabela 1 do documento.
 */
const double tabela1_distanciaDireta[NUM_ESTACOES][NUM_ESTACOES] = {
    /* E1    E2    E3    E4    E5    E6    E7    E8    E9   E10   E11   E12   E13   E14 */
    { 0.0, 10.0, 18.5, 24.8, 36.4, 38.8, 35.8, 25.4, 17.6,  9.1, 16.7, 27.3, 27.6, 29.8}, /* E1  */
    {10.0,  0.0,  8.5, 14.8, 26.6, 29.1, 26.1, 17.3, 10.0,  3.5, 15.5, 20.9, 19.1, 21.8}, /* E2  */
    {18.5,  8.5,  0.0,  6.3, 18.2, 20.6, 17.6, 13.6,  9.4, 10.3, 19.5, 19.1, 12.1, 16.6}, /* E3  */
    {24.8, 14.8,  6.3,  0.0, 12.0, 14.4, 11.5, 12.4, 12.6, 16.7, 23.6, 18.6, 10.6, 15.4}, /* E4  */
    {36.4, 26.6, 18.2, 12.0,  0.0,  3.0,  2.4, 19.4, 23.3, 28.2, 34.2, 24.8, 14.5, 17.9}, /* E5  */
    {38.8, 29.1, 20.6, 14.4,  3.0,  0.0,  3.3, 22.3, 25.7, 30.3, 36.7, 27.6, 15.2, 18.2}, /* E6  */
    {35.8, 26.1, 17.6, 11.5,  2.4,  3.3,  0.0, 20.0, 23.0, 27.3, 34.2, 25.7, 12.4, 15.6}, /* E7  */
    {25.4, 17.3, 13.6, 12.4, 19.4, 22.3, 20.0,  0.0,  8.2, 20.3, 16.1,  6.4, 22.7, 27.6}, /* E8  */
    {17.6, 10.0,  9.4, 12.6, 23.3, 25.7, 23.0,  8.2,  0.0, 13.5, 11.2, 10.9, 21.2, 26.6}, /* E9  */
    { 9.1,  3.5, 10.3, 16.7, 28.2, 30.3, 27.3, 20.3, 13.5,  0.0, 17.6, 24.2, 18.7, 21.2}, /* E10 */
    {16.7, 15.5, 19.5, 23.6, 34.2, 36.7, 34.2, 16.1, 11.2, 17.6,  0.0, 14.2, 31.5, 35.5}, /* E11 */
    {27.3, 20.9, 19.1, 18.6, 24.8, 27.6, 25.7,  6.4, 10.9, 24.2, 14.2,  0.0, 28.8, 33.6}, /* E12 */
    {27.6, 19.1, 12.1, 10.6, 14.5, 15.2, 12.4, 22.7, 21.2, 18.7, 31.5, 28.8,  0.0,  5.1}, /* E13 */
    {29.8, 21.8, 16.6, 15.4, 17.9, 18.2, 15.6, 27.6, 26.6, 21.2, 35.5, 33.6,  5.1,  0.0}  /* E14 */
};

/*
 * TABELA 2: DISTÂNCIAS REAIS PELOS TRILHOS (em km)
 * Utilizada para calcular o custo real de deslocamento g(n).
 * Se o valor for 0.0, não existe conexão ferroviária direta entre as estações.
 * Valores espelhados simetricamente a partir da Tabela 2 do documento.
 */
const double tabela2_distanciaReal[NUM_ESTACOES][NUM_ESTACOES] = {
    /* E1    E2    E3    E4    E5    E6    E7    E8    E9   E10   E11   E12   E13   E14 */
    { 0.0, 10.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, /* E1  */
    {10.0,  0.0,  8.5,  0.0,  0.0,  0.0,  0.0,  0.0, 10.0,  3.5,  0.0,  0.0,  0.0,  0.0}, /* E2  */
    { 0.0,  8.5,  0.0,  6.3,  0.0,  0.0,  0.0,  0.0,  9.4,  0.0,  0.0,  0.0, 18.7,  0.0}, /* E3  */
    { 0.0,  0.0,  6.3,  0.0, 13.0,  0.0,  0.0, 15.3,  0.0,  0.0,  0.0,  0.0, 12.8,  0.0}, /* E4  */
    { 0.0,  0.0,  0.0, 13.0,  0.0,  3.0,  2.4, 30.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, /* E5  */
    { 0.0,  0.0,  0.0,  0.0,  3.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, /* E6  */
    { 0.0,  0.0,  0.0,  0.0,  2.4,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, /* E7  */
    { 0.0,  0.0,  0.0, 15.3, 30.0,  0.0,  0.0,  0.0,  9.6,  0.0,  0.0,  6.4,  0.0,  0.0}, /* E8  */
    { 0.0, 10.0,  9.4,  0.0,  0.0,  0.0,  0.0,  9.6,  0.0,  0.0, 12.2,  0.0,  0.0,  0.0}, /* E9  */
    { 0.0,  3.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, /* E10 */
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 12.2,  0.0,  0.0,  0.0,  0.0,  0.0}, /* E11 */
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  6.4,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, /* E12 */
    { 0.0,  0.0, 18.7, 12.8,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  5.1}, /* E13 */
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  5.1,  0.0}  /* E14 */
};

/*
 * TABELA DE LINHAS: IDENTIFICAÇÃO DAS CORES DAS CONEXÕES
 * 0 = Sem conexão
 * 1 = Linha Azul
 * 2 = Linha Amarela
 * 3 = Linha Vermelha
 * 4 = Linha Verde
 */
const int tabelaLinhas[NUM_ESTACOES][NUM_ESTACOES] = {
    {0,1,0,0,0,0,0,0,0,0,0,0,0,0}, /* E1  */
    {1,0,1,0,0,0,0,0,2,2,0,0,0,0}, /* E2  */
    {0,1,0,1,0,0,0,0,3,0,0,0,3,0}, /* E3  */
    {0,0,1,0,1,0,0,4,0,0,0,0,4,0}, /* E4  */
    {0,0,0,1,0,1,2,2,0,0,0,0,0,0}, /* E5  */
    {0,0,0,0,1,0,0,0,0,0,0,0,0,0}, /* E6  */
    {0,0,0,0,2,0,0,0,0,0,0,0,0,0}, /* E7  */
    {0,0,0,4,2,0,0,0,2,0,0,4,0,0}, /* E8  */
    {0,2,3,0,0,0,0,2,0,0,3,0,0,0}, /* E9  */
    {0,2,0,0,0,0,0,0,0,0,0,0,0,0}, /* E10 */
    {0,0,0,0,0,0,0,0,3,0,0,0,0,0}, /* E11 */
    {0,0,0,0,0,0,0,4,0,0,0,0,0,0}, /* E12 */
    {0,0,3,4,0,0,0,0,0,0,0,0,0,4}, /* E13 */
    {0,0,0,0,0,0,0,0,0,0,0,0,4,0}  /* E14 */
};

/*
 * ESTRUTURA: EstadoBusca
 * Modela os dados de exploração do algoritmo A*.
 * Como a baldeação depende da linha anterior, cada estado no grafo é
 * unicamente representado pelo par: (estação atual, linha do trem).
 */
typedef struct {
    double g[NUM_ESTACOES][NUM_ESTADOS_LINHA];          /* Custo acumulado real da origem até o estado (minutos) */
    double h[NUM_ESTACOES][NUM_ESTADOS_LINHA];          /* Estimativa heurística até o destino (minutos)         */
    double f[NUM_ESTACOES][NUM_ESTADOS_LINHA];          /* Custo total avaliado: f = g + h                      */
    int aberto[NUM_ESTACOES][NUM_ESTADOS_LINHA];        /* 1 se o estado está na fronteira (lista aberta)        */
    int fechado[NUM_ESTACOES][NUM_ESTADOS_LINHA];       /* 1 se o estado já foi explorado (lista fechada)        */
    int paiEstacao[NUM_ESTACOES][NUM_ESTADOS_LINHA];    /* Estação antecessora no melhor caminho                 */
    int paiLinha[NUM_ESTACOES][NUM_ESTADOS_LINHA];      /* Linha utilizada no estado antecessor                  */
    int estacaoDestinoFinal;                            /* Estação de destino alcançada                          */
    int linhaDestinoFinal;                              /* Linha com a qual o destino foi alcançado              */
    double tempoTotalFinal;                             /* Tempo total mínimo em minutos da melhor rota          */
} EstadoBusca;

/*******************************************************************************
 * =============================================================================
 * BLOCO 2: CÁLCULOS DE TEMPO, HEURÍSTICA E TRANSIÇÕES COM BALDEAÇÃO
 * =============================================================================
 * Este bloco implementa as funções responsáveis por traduzir distâncias em
 * minutos, calcular a heurística admissível h(n) a partir da Tabela 1 e
 * computar o custo real das arestas da Tabela 2, incluindo a penalidade de
 * 4 minutos por baldeação quando ocorre troca de linha.
 *******************************************************************************/

/*
 * Converte uma distância em quilômetros para o tempo de deslocamento em minutos.
 *
 * Fórmula física:
 * Velocidade média = 30 km/h = 30 km / 60 min = 0.5 km/min.
 * Tempo (minutos) = Distância (km) / Velocidade (km/min)
 * Tempo = Distância / 0.5 = Distância * 2.0 minutos.
 */
double converterDistanciaParaTempo(double distanciaKm) {
    return (distanciaKm * 60.0) / VELOCIDADE_KM_H;
}

/*
 * Retorna o nome amigável da linha a partir de seu código numérico.
 */
const char* obterNomeLinha(int codigoLinha) {
    switch (codigoLinha) {
        case LINHA_AZUL:     return "Azul";
        case LINHA_AMARELA:  return "Amarela";
        case LINHA_VERMELHA: return "Vermelha";
        case LINHA_VERDE:    return "Verde";
        default:             return "Nenhuma (Partida)";
    }
}

/*
 * Calcula a função heurística h(n):
 * Estimativa do tempo em linha reta entre a estação atual e o destino final.
 * Utiliza os dados da Tabela 1 (distância direta).
 *
 * Propriedade teórica:
 * Como a distância em linha reta é a menor distância geométrica possível entre
 * dois pontos (desigualdade triangular), h(n) nunca superestima o custo real de viagem.
 * Portanto, a heurística é ADMISSÍVEL e CONSISTENTE, garantindo a otimalidade do A*.
 */
double calcularHeuristica(int estacaoAtual, int estacaoDestino) {
    double distanciaLinhaReta = tabela1_distanciaDireta[estacaoAtual][estacaoDestino];
    return converterDistanciaParaTempo(distanciaLinhaReta);
}

/*
 * Calcula o custo real g de transição entre duas estações conectadas.
 *
 * Entradas:
 * - estacaoOrigem: estação de partida do trecho
 * - estacaoDestino: estação de chegada do trecho
 * - linhaAtual: linha na qual o passageiro se encontra antes de embarcar
 * - baldeacaoOcorreu: ponteiro para receber 1 se houve troca de linha, ou 0 caso contrário
 *
 * Regra de baldeação:
 * Se o passageiro já está em uma linha de trem (linhaAtual != 0) e a linha da próxima
 * conexão for diferente (linhaAtual != proximaLinha), adiciona-se o tempo de baldeação
 * de 4.0 minutos.
 */
double calcularCustoAresta(int estacaoOrigem, int estacaoDestino, int linhaAtual, int *baldeacaoOcorreu) {
    int proximaLinha = tabelaLinhas[estacaoOrigem][estacaoDestino];
    double distanciaKm = tabela2_distanciaReal[estacaoOrigem][estacaoDestino];
    double tempoTrecho = converterDistanciaParaTempo(distanciaKm);

    *baldeacaoOcorreu = 0;

    /* Troca de linha em estação intermediária */
    if (linhaAtual != LINHA_NENHUMA && linhaAtual != proximaLinha) {
        tempoTrecho += TEMPO_BALDEACAO_MIN;
        *baldeacaoOcorreu = 1;
    }

    return tempoTrecho;
}

/*******************************************************************************
 * =============================================================================
 * BLOCO 3: ALGORITMO DE BUSCA INFORMADA A* E DIJKSTRA
 * =============================================================================
 * Este bloco implementa o motor de busca do A*.
 * O algoritmo mantém os conjuntos Aberto (fronteira) e Fechado (explorados).
 * A cada iteração, escolhe o nó com menor f(n) = g(n) + h(n), expande seus
 * vizinhos e atualiza os custos se um caminho mais rápido for encontrado.
 * Também aceita o parâmetro usarHeuristica (se for 0, comporta-se como Dijkstra).
 *******************************************************************************/

/*
 * Percorre a fronteira (conjunto aberto) e seleciona o estado (estação, linha)
 * que apresenta o menor valor da função de avaliação f = g + h.
 * Retorna 1 se encontrou um estado, ou 0 se a fronteira estiver vazia.
 */
int selecionarMenorF(const EstadoBusca *busca, int *melhorEstacao, int *melhorLinha) {
    double menorF = VALOR_INFINITO;
    int encontrou = 0;
    int e, l;

    for (e = 0; e < NUM_ESTACOES; e++) {
        for (l = 0; l < NUM_ESTADOS_LINHA; l++) {
            if (busca->aberto[e][l] == 1) {
                if (busca->f[e][l] < menorF) {
                    menorF = busca->f[e][l];
                    *melhorEstacao = e;
                    *melhorLinha = l;
                    encontrou = 1;
                }
            }
        }
    }

    return encontrou;
}

/*
 * Executa o Algoritmo A* (ou Dijkstra se usarHeuristica == 0).
 *
 * Parâmetros:
 * - origem: índice da estação de partida (0 a 13)
 * - destino: índice da estação de chegada (0 a 13)
 * - linhaInicial: linha onde o passageiro começa (0 se ainda não embarcou em nenhuma)
 * - usarHeuristica: 1 para A* (busca informada), 0 para Dijkstra (busca cega)
 * - exibirPassos: 1 para imprimir o log passo a passo da árvore de busca
 * - busca: ponteiro para a estrutura onde os resultados serão salvos
 *
 * Retorna 1 se encontrou caminho, ou 0 caso contrário.
 */
int executarBusca(int origem, int destino, int linhaInicial, int usarHeuristica, int exibirPassos, EstadoBusca *busca) {
    int e, l;
    int estacaoAtual, linhaAtual;
    int vizinho, proximaLinha;
    int passo = 0;

    /* ETAPA 1: Inicialização de todos os estados com valores infinitos */
    for (e = 0; e < NUM_ESTACOES; e++) {
        for (l = 0; l < NUM_ESTADOS_LINHA; l++) {
            busca->g[e][l] = VALOR_INFINITO;
            busca->h[e][l] = 0.0;
            busca->f[e][l] = VALOR_INFINITO;
            busca->aberto[e][l] = 0;
            busca->fechado[e][l] = 0;
            busca->paiEstacao[e][l] = -1;
            busca->paiLinha[e][l] = -1;
        }
    }

    busca->estacaoDestinoFinal = -1;
    busca->linhaDestinoFinal = -1;
    busca->tempoTotalFinal = VALOR_INFINITO;

    /* ETAPA 2: Configuração e inserção do nó raiz (origem) na lista aberta */
    busca->g[origem][linhaInicial] = 0.0;
    busca->h[origem][linhaInicial] = usarHeuristica ? calcularHeuristica(origem, destino) : 0.0;
    busca->f[origem][linhaInicial] = busca->g[origem][linhaInicial] + busca->h[origem][linhaInicial];
    busca->aberto[origem][linhaInicial] = 1;

    if (exibirPassos) {
        printf("\n=======================================================================\n");
        printf(" EXECUTANDO BUSCA %s (Origem: E%d -> Destino: E%d)\n",
               usarHeuristica ? "A* (INFORMADA)" : "DIJKSTRA (CEGA)",
               origem + 1, destino + 1);
        printf("=======================================================================\n");
    }

    /* ETAPA 3: Laço principal de busca */
    while (selecionarMenorF(busca, &estacaoAtual, &linhaAtual)) {
        passo++;

        /* Retira o estado atual da lista aberta e transfere para a fechada */
        busca->aberto[estacaoAtual][linhaAtual] = 0;
        busca->fechado[estacaoAtual][linhaAtual] = 1;

        if (exibirPassos) {
            printf("Passo %2d: Expandindo E%-2d | Linha: %-8s | g = %5.1f min | h = %5.1f min | f = %5.1f min\n",
                   passo, estacaoAtual + 1, obterNomeLinha(linhaAtual),
                   busca->g[estacaoAtual][linhaAtual],
                   busca->h[estacaoAtual][linhaAtual],
                   busca->f[estacaoAtual][linhaAtual]);
        }

        /* ETAPA 4: Verificação do teste de meta (chegada ao destino) */
        if (estacaoAtual == destino) {
            busca->estacaoDestinoFinal = estacaoAtual;
            busca->linhaDestinoFinal = linhaAtual;
            busca->tempoTotalFinal = busca->g[estacaoAtual][linhaAtual];

            if (exibirPassos) {
                printf("-----------------------------------------------------------------------\n");
                printf("-> Meta alcancada com sucesso na Estacao E%d! Tempo total: %.2f minutos\n",
                       destino + 1, busca->tempoTotalFinal);
                printf("=======================================================================\n");
            }
            return 1;
        }

        /* ETAPA 5: Expansão dos vizinhos conectados na Tabela 2 */
        for (vizinho = 0; vizinho < NUM_ESTACOES; vizinho++) {
            /* Verifica se há trilhos conectando as estações */
            if (tabela2_distanciaReal[estacaoAtual][vizinho] > 0.0) {
                int baldeacao = 0;
                double custoAresta;
                double novoG, novoH, novoF;

                proximaLinha = tabelaLinhas[estacaoAtual][vizinho];
                custoAresta = calcularCustoAresta(estacaoAtual, vizinho, linhaAtual, &baldeacao);
                novoG = busca->g[estacaoAtual][linhaAtual] + custoAresta;

                /* Se encontrou um caminho estritamente mais rápido para este estado vizinho */
                if (novoG < busca->g[vizinho][proximaLinha]) {
                    novoH = usarHeuristica ? calcularHeuristica(vizinho, destino) : 0.0;
                    novoF = novoG + novoH;

                    busca->g[vizinho][proximaLinha] = novoG;
                    busca->h[vizinho][proximaLinha] = novoH;
                    busca->f[vizinho][proximaLinha] = novoF;
                    busca->paiEstacao[vizinho][proximaLinha] = estacaoAtual;
                    busca->paiLinha[vizinho][proximaLinha] = linhaAtual;

                    /* Reabre o nó se já estivesse fechado, ou insere na lista aberta */
                    busca->fechado[vizinho][proximaLinha] = 0;
                    busca->aberto[vizinho][proximaLinha] = 1;

                    if (exibirPassos) {
                        printf("         -> Gerou vizinho E%-2d (Linha %-8s) | Custo aresta: %4.1f min%s | Novo f = %5.1f\n",
                               vizinho + 1, obterNomeLinha(proximaLinha), custoAresta,
                               baldeacao ? " [Troca de Linha +4m]" : "", novoF);
                    }
                }
            }
        }
    }

    return 0; /* Nenhum caminho encontrado */
}

/*******************************************************************************
 * =============================================================================
 * BLOCO 4: RECONSTRUÇÃO DO PERCURSO E RELATÓRIO DETALHADO DA ROTA
 * =============================================================================
 * Este bloco é responsável por navegar de volta pelos nós pais desde o destino
 * até a origem, organizando a ordem cronológica da viagem e imprimindo uma
 * tabela formatada e clara com as estações, linhas, baldeações e tempos.
 *******************************************************************************/

/*
 * Reconstrói o trajeto voltando do destino até a origem através dos ponteiros
 * de predecessores (paiEstacao e paiLinha).
 * Retorna a quantidade de paradas no percurso.
 */
int reconstruirCaminho(const EstadoBusca *busca, int origem, int destino,
                       int estacoesCaminho[], int linhasCaminho[]) {
    int estacaoAtual = destino;
    int linhaAtual = busca->linhaDestinoFinal;
    int temporarioEstacoes[NUM_ESTACOES * NUM_ESTADOS_LINHA];
    int temporarioLinhas[NUM_ESTACOES * NUM_ESTADOS_LINHA];
    int totalPassos = 0;
    int i;

    if (busca->tempoTotalFinal >= VALOR_INFINITO || linhaAtual == -1) {
        return 0;
    }

    /* Retrocede do destino até a estação de origem */
    while (estacaoAtual != -1 && linhaAtual != -1) {
        temporarioEstacoes[totalPassos] = estacaoAtual;
        temporarioLinhas[totalPassos] = linhaAtual;
        totalPassos++;

        if (estacaoAtual == origem) {
            break;
        }

        int paiE = busca->paiEstacao[estacaoAtual][linhaAtual];
        int paiL = busca->paiLinha[estacaoAtual][linhaAtual];

        estacaoAtual = paiE;
        linhaAtual = paiL;
    }

    /* Inverte os vetores para obter a ordem cronológica da partida à chegada */
    for (i = 0; i < totalPassos; i++) {
        estacoesCaminho[i] = temporarioEstacoes[totalPassos - 1 - i];
        linhasCaminho[i] = temporarioLinhas[totalPassos - 1 - i];
    }

    return totalPassos;
}

/*
 * Apresenta na tela o relatório completo, limpo e legível do percurso ótimo.
 */
void exibirRelatorioPercurso(const EstadoBusca *busca, int origem, int destino) {
    int estacoesCaminho[NUM_ESTACOES * NUM_ESTADOS_LINHA];
    int linhasCaminho[NUM_ESTACOES * NUM_ESTADOS_LINHA];
    int totalPassos;
    int i;
    double distanciaTotalKm = 0.0;
    int totalBaldeacoes = 0;
    double tempoAcumulado = 0.0;

    totalPassos = reconstruirCaminho(busca, origem, destino, estacoesCaminho, linhasCaminho);

    if (totalPassos == 0) {
        printf("\n[AVISO] Nao foi possivel encontrar um trajeto entre E%d e E%d.\n", origem + 1, destino + 1);
        return;
    }

    printf("\n===============================================================================\n");
    printf("                  RELATORIO DO MELHOR TRAJETO DE METRO\n");
    printf("===============================================================================\n");
    printf(" Partida : Estacao E%-2d\n", origem + 1);
    printf(" Destino : Estacao E%-2d\n", destino + 1);
    printf("-------------------------------------------------------------------------------\n");
    printf(" #  | Trecho          | Linha Usada     | Distancia | Tempo Trecho | Tempo Acum.\n");
    printf("-------------------------------------------------------------------------------\n");

    for (i = 0; i < totalPassos; i++) {
        int eAtual = estacoesCaminho[i];

        if (i == 0) {
            printf("%2d  | Embarque em E%-2d | %-15s | %6.1f km |   %5.1f min  |  %5.1f min\n",
                   i + 1, eAtual + 1, "Ponto de Partida", 0.0, 0.0, 0.0);
        } else {
            int eAnterior = estacoesCaminho[i - 1];
            int linhaTrecho = linhasCaminho[i];
            int linhaAnterior = linhasCaminho[i - 1];
            double distKm = tabela2_distanciaReal[eAnterior][eAtual];
            double tempoDeslocamento = converterDistanciaParaTempo(distKm);
            int houveBaldeacao = (linhaAnterior != LINHA_NENHUMA && linhaAnterior != linhaTrecho);

            distanciaTotalKm += distKm;

            if (houveBaldeacao) {
                totalBaldeacoes++;
                tempoAcumulado += TEMPO_BALDEACAO_MIN;
                printf("    | [BALDEACAO]     | Troca p/ %-7s |     -     |  + %4.1f min  |  %5.1f min\n",
                       obterNomeLinha(linhaTrecho), TEMPO_BALDEACAO_MIN, tempoAcumulado);
            }

            tempoAcumulado += tempoDeslocamento;

            printf("%2d  | E%-2d -> E%-2d       | Linha %-9s | %6.1f km |   %5.1f min  |  %5.1f min\n",
                   i + 1, eAnterior + 1, eAtual + 1, obterNomeLinha(linhaTrecho),
                   distKm, tempoDeslocamento, tempoAcumulado);
        }
    }

    printf("===============================================================================\n");
    printf(" RESUMO DA VIAGEM:\n");
    printf(" - Distancia total percorrida nos trilhos: %.1f km\n", distanciaTotalKm);
    printf(" - Quantidade de baldeacoes necessarias : %d (4 min cada)\n", totalBaldeacoes);
    printf(" - Tempo total estimado de viagem       : %.2f minutos (%.1f h)\n",
           tempoAcumulado, tempoAcumulado / 60.0);
    printf("===============================================================================\n");
}

/*******************************************************************************
 * =============================================================================
 * BLOCO 5: INTERFACE COM USUÁRIO, DEMONSTRAÇÃO E QUESTÕES TEÓRICAS (MAIN)
 * =============================================================================
 * Este bloco é responsável pelo menu interativo de opções, execução do caso
 * solicitado no documento (E6 Linha Azul para E13 Linha Vermelha), exibição das
 * respostas teóricas exigidas e controle do fluxo principal da aplicação.
 *******************************************************************************/

/*
 * Exibe um mapa esquemático em texto resumindo as 14 estações e suas linhas.
 */
void exibirMapaEsquematico(void) {
    printf("\n=======================================================================\n");
    printf("              MAPA ESQUEMATICO DAS LINHAS DO METRO DE PARIS\n");
    printf("=======================================================================\n");
    printf(" Linha 1 (Azul)    : E1 <-> E2 <-> E3 <-> E4 <-> E5 <-> E6\n");
    printf(" Linha 2 (Amarela) : E10 <-> E2 <-> E9 <-> E8 <-> E5 <-> E7\n");
    printf(" Linha 3 (Vermelha): E11 <-> E9 <-> E3 <-> E13\n");
    printf(" Linha 4 (Verde)   : E12 <-> E8 <-> E4 <-> E13 <-> E14\n");
    printf("-----------------------------------------------------------------------\n");
    printf(" Estacoes de Integracao (Baldeacao):\n");
    printf(" - E2 : Linha 1 (Azul) e Linha 2 (Amarela)\n");
    printf(" - E3 : Linha 1 (Azul) e Linha 3 (Vermelha)\n");
    printf(" - E4 : Linha 1 (Azul) e Linha 4 (Verde)\n");
    printf(" - E5 : Linha 1 (Azul) e Linha 2 (Amarela)\n");
    printf(" - E8 : Linha 2 (Amarela) e Linha 4 (Verde)\n");
    printf(" - E9 : Linha 2 (Amarela) e Linha 3 (Vermelha)\n");
    printf(" - E13: Linha 3 (Vermelha) e Linha 4 (Verde)\n");
    printf("=======================================================================\n");
}

/*
 * Exibe na tela as respostas completas das duas questões solicitadas
 * no arquivo ProjetoGrafos.doc.
 */
void exibirRespostasTeoricas(void) {
    printf("\n=======================================================================\n");
    printf("      RESPOSTAS DAS QUESTOES DO TRABALHO (ProjetoGrafos.doc)\n");
    printf("=======================================================================\n\n");

    printf("--- QUESTAO 1: Descricao do funcionamento do algoritmo em linguagem corrente ---\n");
    printf("1. O algoritmo A* calcula para cada no a funcao de avaliacao f(n) = g(n) + h(n),\n");
    printf("2. onde g(n) e o tempo real acumulado desde a origem e h(n) e a estimativa em\n");
    printf("3. linha reta ate o destino, garantindo uma busca informada e admissivel.\n");
    printf("4. A cada passo, o algoritmo seleciona na fronteira o estado com menor f(n),\n");
    printf("5. somando 4 minutos por baldeacao ate expandir o destino com o menor tempo total.\n\n");

    printf("--- QUESTAO 2: Arvore de busca desenvolvida pelo A* (Ei = E6 Azul -> Ef = E13) ---\n");
    printf("1. Raiz: Expandir E6 (Linha Azul, g=0.0, h=30.4, f=30.4)\n");
    printf("   -> Gera vizinho E5 (Linha Azul, g=6.0, h=29.0, f=35.0)\n");
    printf("2. Expandir E5 (Linha Azul, f=35.0)\n");
    printf("   -> Gera E7 (Linha Amarela, troca de linha: g=6+4+4.8=14.8, h=24.8, f=39.6)\n");
    printf("   -> Gera E4 (Linha Azul, mesma linha: g=6+26.0=32.0, h=21.2, f=53.2)\n");
    printf("   -> Gera E8 (Linha Amarela, troca de linha: g=6+4+60.0=70.0, h=45.4, f=115.4)\n");
    printf("3. Menor f na fronteira: Expandir E7 (Linha Amarela, f=39.6)\n");
    printf("   -> Estacao terminal sem outras conexoes (beco sem saida).\n");
    printf("4. Menor f na fronteira: Expandir E4 (Linha Azul, f=53.2)\n");
    printf("   -> Gera E13 (Linha Verde, troca de linha: g=32+4+25.6=61.6, h=0.0, f=61.6)\n");
    printf("   -> Gera E3  (Linha Azul, mesma linha: g=32+12.6=44.6, h=24.2, f=68.8)\n");
    printf("   -> Gera E8  (Linha Verde, troca de linha: g=32+4+30.6=66.6, h=45.4, f=112.0)\n");
    printf("5. Menor f na fronteira: Expandir E13 (f=61.6) -> Meta alcancada!\n");
    printf("   -> Percurso otimo encontrado: E6 -> E5 -> E4 [Baldeacao] -> E13\n");
    printf("   -> Tempo Total: 61.60 minutos\n");
    printf("=======================================================================\n");
}

/*
 * Executa a demonstração do caso de teste padrão exigido no documento:
 * Estação inicial: E6 (Linha Azul)
 * Estação final: E13 (Linha Vermelha / Estação 13)
 */
void executarCasoTrabalho(void) {
    EstadoBusca busca;
    int origem = 5;       /* E6 (indice 5) */
    int destino = 12;     /* E13 (indice 12) */
    int linhaInicial = LINHA_AZUL; /* Linha Azul = 1 */

    printf("\n>>> EXECUTANDO CASO DE TESTE ESPECIFICADO NO PROJETO <<<\n");
    printf("Origem : E6 (Linha Azul)\n");
    printf("Destino: E13 (Linha Vermelha)\n");

    if (executarBusca(origem, destino, linhaInicial, 1, 1, &busca)) {
        exibirRelatorioPercurso(&busca, origem, destino);
    } else {
        printf("\nErro ao calcular a rota para o caso de teste.\n");
    }
}

/*
 * Executa uma consulta interativa onde o usuário escolhe qualquer estação
 * de origem e de destino.
 */
void executarConsultaInterativa(void) {
    EstadoBusca busca;
    int origemUsuario, destinoUsuario;
    int linhaInicial = LINHA_NENHUMA;
    int exibirPassos;

    printf("\nDigite a estacao de ORIGEM (1 a 14): ");
    if (scanf("%d", &origemUsuario) != 1 || origemUsuario < 1 || origemUsuario > 14) {
        printf("Entrada invalida! Retornando ao menu principal...\n");
        while (getchar() != '\n');
        return;
    }

    printf("Digite a estacao de DESTINO (1 a 14): ");
    if (scanf("%d", &destinoUsuario) != 1 || destinoUsuario < 1 || destinoUsuario > 14) {
        printf("Entrada invalida! Retornando ao menu principal...\n");
        while (getchar() != '\n');
        return;
    }

    printf("Deseja ver a expansao passo a passo da arvore A*? (1 = Sim, 0 = Nao): ");
    if (scanf("%d", &exibirPassos) != 1) {
        exibirPassos = 0;
        while (getchar() != '\n');
    }

    /* Converte de 1-based (usuário) para 0-based (vetores) */
    int origIdx = origemUsuario - 1;
    int destIdx = destinoUsuario - 1;

    if (executarBusca(origIdx, destIdx, linhaInicial, 1, exibirPassos, &busca)) {
        exibirRelatorioPercurso(&busca, origIdx, destIdx);
    } else {
        printf("\nNao foi possivel encontrar um trajeto entre E%d e E%d.\n", origemUsuario, destinoUsuario);
    }
}

/*
 * Demonstração comparativa entre o A* (Busca Informada) e Dijkstra (Busca Cega).
 * Mostra como a heurística reduz o número de passos explorados.
 */
void compararAEstrelaEDijkstra(void) {
    EstadoBusca buscaAEstrela;
    EstadoBusca buscaDijkstra;
    int origIdx = 5;  /* E6 */
    int destIdx = 12; /* E13 */

    printf("\n=======================================================================\n");
    printf("        COMPARATIVO: BUSCA INFORMADA (A*) vs BUSCA CEGA (DIJKSTRA)\n");
    printf("                    Trajeto: E6 ate E13\n");
    printf("=======================================================================\n");

    printf("\n--- 1. Executando Algoritmo A* (com heuristica h da Tabela 1) ---\n");
    executarBusca(origIdx, destIdx, LINHA_AZUL, 1, 1, &buscaAEstrela);

    printf("\n--- 2. Executando Algoritmo de Dijkstra (sem heuristica, h = 0) ---\n");
    executarBusca(origIdx, destIdx, LINHA_AZUL, 0, 1, &buscaDijkstra);

    printf("\n=======================================================================\n");
    printf(" CONCLUSAO COMPARATIVA:\n");
    printf(" - Ambos os algoritmos encontram a mesma rota otima (Tempo: %.2f min).\n", buscaAEstrela.tempoTotalFinal);
    printf(" - O A* utiliza a estimativa em linha reta para guiar a busca diretamente\n");
    printf("   em direcao ao destino, evitando explorar nos geograficamente distantes.\n");
    printf(" - O algoritmo de Dijkstra expande nos em todas as direcoes igualmente.\n");
    printf("=======================================================================\n");
}

/*
 * Função principal (ponto de entrada do programa).
 */
int main(void) {
    int opcao;

    do {
        printf("\n=======================================================================\n");
        printf("               PAINEL DO METRO DE PARIS - SISTEMA DE ROTAS\n");
        printf("=======================================================================\n");
        printf(" 1 - Executar caso do trabalho: E6 (Azul) -> E13 (Vermelha)\n");
        printf(" 2 - Calcular melhor rota entre quaisquer duas estacoes\n");
        printf(" 3 - Comparar Algoritmo A* (Informada) vs Dijkstra (Cega)\n");
        printf(" 4 - Ver mapa esquematico das estacoes e conexoes\n");
        printf(" 5 - Exibir respostas das questoes teoricas do trabalho\n");
        printf(" 0 - Sair do programa\n");
        printf("-----------------------------------------------------------------------\n");
        printf(" Escolha uma opcao: ");

        if (scanf("%d", &opcao) != 1) {
            printf("\nOpcao invalida! Digite um numero.\n");
            while (getchar() != '\n');
            opcao = -1;
            continue;
        }

        switch (opcao) {
            case 1:
                executarCasoTrabalho();
                break;
            case 2:
                executarConsultaInterativa();
                break;
            case 3:
                compararAEstrelaEDijkstra();
                break;
            case 4:
                exibirMapaEsquematico();
                break;
            case 5:
                exibirRespostasTeoricas();
                break;
            case 0:
                printf("\nEncerrando o sistema do Metro de Paris. Ate logo!\n");
                break;
            default:
                printf("\nOpcao invalida! Escolha uma opcao entre 0 e 5.\n");
                break;
        }

    } while (opcao != 0);

    return 0;
}
