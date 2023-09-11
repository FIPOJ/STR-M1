#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <immintrin.h>
#include <ammintrin.h>
#include <stdint.h>
#include <windows.h>

#include "pgmfiles.h"
#include "diff2d.h"



//gcc -o fda pgmtolist.c pgmfiles.c diff2d.c main.c -lm

// Adição: Definir o número de threads
#define NUM_THREADS 4

#define LUT_SIZE 1000

float LUT[LUT_SIZE];  // Tabela de Pesquisa (LUT)

// Função para inicializar o LUT
void initialize_LUT() {
  for (int i = 0; i < LUT_SIZE; ++i) {
    LUT[i] = pow(i, 2) + sqrt(i) - log(i);
  }
}

// Nova função para o algoritmo FDA com uso de LUT
void diff2d_LUT(float **matrix, int nLinhas, int nColunas) {
    float LUT[256];  // A tabela de pesquisa, ajuste o tamanho conforme necessário

    // Preencher a LUT. Aqui é só um exemplo; use seus próprios cálculos.
    for (int i = 0; i < 256; ++i) {
        LUT[i] = sqrt((float) i); // Substitua esta operação pelo seu cálculo real
    }

    // Aplicar a LUT na matriz
    for (int i = 0; i < nLinhas; ++i) {
        for (int j = 0; j < nColunas; ++j) {
            int valor = (int) matrix[i][j]; // Assumindo que os valores estão entre 0 e 255
            if (valor >= 0 && valor < 256) {
                matrix[i][j] = LUT[valor];
            }
        }
    }
}

DWORD WINAPI MinhaThreadLUT(LPVOID param) {
    // Converter o argumento de LPVOID para o tipo de dados apropriado
    int* thread_data = (int*)param;

    // Obter os dados da matriz a partir do argumento
    float **matrix = (float **) thread_data[0];
    int tamanho_matriz = thread_data[1];
    int inicio = thread_data[2];
    int fim = thread_data[3];
    int nColunas = thread_data[4];

    // Processar a parte da matriz correspondente usando LUT
    for (int i = inicio; i < fim; ++i) {
        for (int j = 0; j < nColunas; ++j) {
            int valor = (int) matrix[i][j];
            if (valor >= 0 && valor < LUT_SIZE) { // assumindo que os valores estão entre 0 e LUT_SIZE - 1
                matrix[i][j] = LUT[valor];
            }
        }
    }

    return 0;
}


DWORD WINAPI MinhaThread(LPVOID param) {
    // Converter o argumento de LPVOID para o tipo de dados apropriado
    int* thread_data = (int*)param;

    // Obter os dados da matriz a partir do argumento
    int* matriz = thread_data[0];
    int tamanho_matriz = thread_data[1];
    int inicio = thread_data[2];
    int fim = thread_data[3];

    // Processar a parte da matriz correspondente
    for (int i = inicio; i < fim; i++) {
        // Realizar o processamento necessário
        matriz[i] *= 2;  // Exemplo simples: multiplicar todos os elementos por 2
    }

    return 0;
}


// Função para ler o contador de ciclos
static inline uint64_t rdtsc() {
    uint32_t lo, hi;
    __asm__ volatile (
        "rdtsc"
        : "=a" (lo), "=d" (hi)
    );
    return ((uint64_t)hi << 32) | lo;
}


void main (int argc, char **argv) {
  char   row[80];
  float  **matrix;
  int i, j;
  FILE   *inimage, *outimage;
  long   imax;
  float  lambda;
  int result;
  eightBitPGMImage *PGMImage;
  time_t begin, end, begin_total, end_total;
  double time_spent_total, time_spent_in_arq, time_spent_out_arq, time_spent_diff2d;
  uint64_t start_cycles, end_cycles, cycles_in_arq, cycles_out_arq, cycles_diff2d, cycle_begin, cycle_end;


  /* ---- read image name  ---- */

  begin_total = clock();
  start_cycles = rdtsc();

  PGMImage = (eightBitPGMImage *) malloc(sizeof(eightBitPGMImage));

  if (!argv[1])
  {
    printf("name of input PGM image file (with extender): ");
    scanf("%s", PGMImage->fileName);
  }
  else
  {
    strcpy(PGMImage->fileName, argv[1]);
  }

  begin = clock();
  cycle_begin = rdtsc();

  result = read8bitPGM(PGMImage);

  if(result < 0)
    {
      printPGMFileError(result);
      exit(result);
    }

  /* ---- allocate storage for matrix ---- */

  matrix = (float **) malloc (PGMImage->x * sizeof(float *));
  if (matrix == NULL)
    {
      printf("not enough storage available\n");
      exit(1);
    }
  for (i=0; i<PGMImage->x; i++)
    {
      matrix[i] = (float *) malloc (PGMImage->y * sizeof(float));
      if (matrix[i] == NULL)
        {
          printf("not enough storage available\n");
          exit(1);
        }
    }

  /* ---- read image data into matrix ---- */

 for (i=0; i<PGMImage->x; i++)
    for (j=0; j<PGMImage->y; j++)
      matrix[i][j] = (float) *(PGMImage->imageData + (i*PGMImage->y) + j);

  end = clock();
  cycle_end = rdtsc();

  time_spent_in_arq = ((double) (end - begin)) / CLOCKS_PER_SEC;
  cycles_in_arq = cycle_end - cycle_begin;


  /* ---- process image ---- */

  printf("contrast paramter lambda (>0) : ");
  //~ gets(row);  sscanf(row, "%f", &lambda);
  scanf("%f", &lambda);
  printf("number of iterations: ");
  //~ gets(row);  sscanf(row, "%ld", &imax);
  scanf("%ld", &imax);

  begin = clock();
  cycle_begin = rdtsc();

  for (i=1; i<=imax; i++)
    {
      printf("iteration number: %3ld \n", i);
      diff2d (0.5, lambda, PGMImage->x, PGMImage->y, matrix);
    }

  end = clock();
  cycle_end = rdtsc();

  time_spent_diff2d = ((double) (end - begin)) / CLOCKS_PER_SEC;
  cycles_diff2d = cycle_end - cycle_begin;

  /* copy the Result Image to PGM Image/File structure */

  for (i=0; i<PGMImage->x; i++)
    for (j=0; j<PGMImage->y; j++)
      *(PGMImage->imageData + i*PGMImage->y + j) = (char) matrix[i][j];

  /* ---- write image ---- */

/*
  if (!argv[2])
  {
    printf("name of output PGM image file (with extender): ");
    scanf("%s", PGMImage->fileName);
  }

  else
  {
    strcpy(PGMImage->fileName, argv[2]);
  }
  */

  printf("name of output PGM image file (with extender): ");
  scanf("%s", PGMImage->fileName);

  begin = clock();
  cycle_begin = rdtsc();

  write8bitPGM(PGMImage);

  /* ---- disallocate storage ---- */

  for (i=0; i<PGMImage->x; i++)
    free(matrix[i]);
  free(matrix);

  free(PGMImage->imageData);
  free(PGMImage);

  end = clock();
  cycle_end = rdtsc();

  end_total = clock();
  end_cycles = rdtsc();


  time_spent_out_arq   = ((double) (end - begin)) / CLOCKS_PER_SEC;
  cycles_out_arq = cycle_end - cycle_begin;
  time_spent_total   = ((double) (end_total - begin_total)) / CLOCKS_PER_SEC;


  printf("\nTempo de execucao usando bicliotecas proprias do windows:\n");
  printf("\nTempo de CPU usado na abertura do arquivo: %lf segundos", time_spent_in_arq);
  printf("\nTempo de CPU usado na execucao do FDA: %lf segundos", time_spent_diff2d);
  printf("\nTempo de CPU usado no fechamento do arquivo: %lf segundos", time_spent_out_arq);
  printf("\nTempo de CPU usado no total da execucao: %lf segundos", time_spent_total);


  printf("\n\nTempo de execucao usando contagem de ciclos:\n");
  printf("\nNumero de ciclos de CPU na abertura do arquivo: %llu\n", cycles_in_arq);
  printf("\nNumero de ciclos de CPU na execucao do FDA: %llu\n", cycles_diff2d);
  printf("\nNumero de ciclos de CPU no fechamento do arquivo: %llu\n", cycles_out_arq);
  printf("\nNumero de ciclos de CPU no total: %llu\n", end_cycles - start_cycles);

  HANDLE threads[NUM_THREADS];
  DWORD threadID;
  int step = PGMImage->x / NUM_THREADS;

  // Para armazenar tempos e ciclos
  double time_spent_parallel[NUM_THREADS];
  uint64_t cycles_parallel[NUM_THREADS];
  double total_time_spent_parallel = 0.0;
  uint64_t total_cycles_parallel = 0;
  printf("\n| Iteration |     Tempo (s)     | Ciclos CPU |\n");
  printf("|-----------|-------------------|------------|\n");

  // Processamento de imagem (parte original)
  for (i=1; i<=imax; i++) {
    uint64_t start_cycles_parallel, end_cycles_parallel;  // Para medir os ciclos da seção paralela
    time_t begin_parallel, end_parallel;
    printf("| %9ld | ", i);  // Mostrar o número da iteração na tabela

    // Início da medição de tempo e ciclos para a seção paralela
    begin_parallel = clock();
    start_cycles_parallel = rdtsc();

    // Adição: Criar threads para paralelizar o processamento
    for (int j = 0; j < NUM_THREADS; j++) {
      int* thread_data = malloc(4 * sizeof(int));
      thread_data[0] = (int) matrix;
      thread_data[1] = PGMImage->x;
      thread_data[2] = j * step;
      thread_data[3] = (j + 1) * step;
      threads[j] = CreateThread(NULL, 0, MinhaThread, thread_data, 0, &threadID);
    }

    // Adição: Aguardar todas as threads terminarem
    for (int j = 0; j < NUM_THREADS; j++) {
      WaitForSingleObject(threads[j], INFINITE);
      CloseHandle(threads[j]);
    }

    // Fim da medição de tempo e ciclos para a seção paralela
    end_parallel = clock();
    end_cycles_parallel = rdtsc();

    time_spent_parallel[i-1] = ((double) (end_parallel - begin_parallel)) / CLOCKS_PER_SEC;
    cycles_parallel[i-1] = end_cycles_parallel - start_cycles_parallel;
    total_time_spent_parallel += time_spent_parallel[i-1];
    total_cycles_parallel += cycles_parallel[i-1];
    printf("%17.10lf | %10llu |\n", time_spent_parallel[i-1], cycles_parallel[i-1]);
  }
  printf("Tempo total: %17.10lf segundos\n", total_time_spent_parallel);
  printf("Ciclos totais: %llu\n", total_cycles_parallel);

  initialize_LUT(); // Inicializa o LUT

  // Loop para as iterações do FDA com LUT
  imax = 10; // ou qualquer número de iterações que você queira
  clock_t start_time_lut = clock();
  uint64_t   cycle_begin_lut = rdtsc();
  for (int i = 1; i <= imax; i++) {
    printf("iteration number with LUT: %3d \n", i);
    diff2d_LUT(matrix, PGMImage->x, PGMImage->y);
  }
  clock_t end_time_lut = clock();
  uint64_t cycle_end_lut = rdtsc();

  double time_spent_lut = ((double) (end_time_lut - start_time_lut)) / CLOCKS_PER_SEC;
  uint64_t cycles_lut = cycle_end_lut - cycle_begin_lut;

  // Estimativa de ciclos de CPU
  printf("\nTempo de CPU usado no total da execucao: %lf segundos", time_spent_lut);
  printf("\nNumero de ciclos de CPU no total: %llu\n", cycle_end_lut - cycle_begin_lut);

  HANDLE threadsLUT[NUM_THREADS];
  DWORD threadIDLUT;
  step = PGMImage->x / NUM_THREADS;

  // Inicialize o LUT aqui, se ainda não tiver feito
  initialize_LUT();

  // Início da medição de tempo e ciclos para a seção paralela com LUT
  clock_t begin_lut_parallel = clock();
  uint64_t start_cycles_lut_parallel = rdtsc();

  // Criar threads para paralelizar o processamento com LUT
  for (int j = 0; j < NUM_THREADS; j++) {
    int* thread_data = malloc(5 * sizeof(int));
    thread_data[0] = (int)matrix;
    thread_data[1] = PGMImage->x;
    thread_data[2] = j * step;
    thread_data[3] = (j + 1) * step;
    thread_data[4] = PGMImage->y;
    threadsLUT[j] = CreateThread(NULL, 0, MinhaThreadLUT, thread_data, 0, &threadIDLUT);
  }

  // Aguardar todas as threads terminarem
  for (int j = 0; j < NUM_THREADS; j++) {
    WaitForSingleObject(threadsLUT[j], INFINITE);
    CloseHandle(threadsLUT[j]);
  }

  // Fim da medição de tempo e ciclos para a seção paralela com LUT
  clock_t end_lut_parallel = clock();
  uint64_t end_cycles_lut_parallel = rdtsc();

  double time_spent_lut_parallel = ((double) (end_lut_parallel - begin_lut_parallel)) / CLOCKS_PER_SEC;
  uint64_t cycles_lut_parallel = end_cycles_lut_parallel - start_cycles_lut_parallel;

  // Mostrar os resultados
  printf("\nTempo de CPU usando paralelismo e LUT: %lf segundos", time_spent_lut_parallel);
  printf("\nNúmero de ciclos de CPU usando paralelismo e LUT: %llu\n", cycles_lut_parallel);
}
