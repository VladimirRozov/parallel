#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>

const int A = 321;
int *percent;
pthread_mutex_t print_mutex, chunk_mutex, percent_mutex;
double *m1, *m2, *m2_copy;
int N, FOR_I, THREAD_NUM, start_chunk_size;

//Барьер - глобальная переменная
static pthread_barrier_t barrier;

struct main_args {
    int index;
    int id;
};

struct timing {
    double start;
    double gen;
    double map;
    double merge;
    double sort;
    double reduce;
    double end;
};

void print_timing(struct timing *timing) {
    printf("Start time: %f ms\n", timing->start);
    printf("Gen time: %f ms\n", timing->gen - timing->start);
    //printf("Generation time: %f sec\n", timing->gen);
    printf("Map time: %f ms\n", timing->map - timing->gen);
    //printf("Map time: %f sec\n", timing->map);
    printf("Merge time: %f ms\n", timing->merge - timing->map);
   // printf("Merge time: %f sec\n", timing->merge);
    printf("Sort time: %f ms\n", timing->sort - timing->merge);
   // printf("Sort time: %f sec\n", timing->sort);
    //printf("Reduce time: %f sec\n", timing->reduce);
    printf("Reduce time: %f ms\n", timing->reduce - timing->sort);
    printf("End time: %f ms\n", timing->end);
    printf("All time: %f ms\n", timing->end - timing->start);
}

double get_wtime() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + ((double)t.tv_usec) / 1000.0;
}

int cth_index = 0;
int cbrtpi_index = 0;
int merge_index = 0;

//гиперболический котангенс корня числа
double cth(int x){
    return (exp(2*sqrt(x))+1)/(exp(2*sqrt(x))-1);
}

//Кубический корень после умножения на число Пи
double cbrtpi(int x){
    return cbrt(x * M_PI);
}

int get_cth_index() {
    int r;
    pthread_mutex_lock(&chunk_mutex);
    r = cth_index;
    cth_index += start_chunk_size;
    if (r >= N) {
        r = -1; // stop
    }
    pthread_mutex_unlock(&chunk_mutex);
    return r;
}

int get_cbrtpi_index() {
    int r;
    pthread_mutex_lock(&chunk_mutex);
    r = cbrtpi_index;
    cbrtpi_index += start_chunk_size;
    if (r >= N / 2) {
        r = -1; // stop
    }
    pthread_mutex_unlock(&chunk_mutex);
    return r;
}

int get_merge_index() {
    int r;
    pthread_mutex_lock(&chunk_mutex);
    r = merge_index;
    merge_index += start_chunk_size;
    if (r >= N / 2) {
        r = -1; // stop
    }
    pthread_mutex_unlock(&chunk_mutex);
    return r;
}

//глупая сортировка
double stupid_sort() {
   int i;
   for(i = 0; i < N/2; i++) {
      if(m2[i] > m2[i + 1]){
         int tmp = m2[i];
         m2[i] = m2[i + 1];
         m2[i + 1] = tmp;
         i = -1;
      }
   }
   return 0;
}

void cth_part(int max_size) {
    int start_i = get_cth_index();

    while (start_i >= 0) {
        int len = start_chunk_size;
        if (start_i + len > max_size) len = max_size - start_i;

        int counter = 0;
        while (counter < len) {
            m1[start_i + counter] = cth(m1[start_i + counter]);
            ++counter;
        }
        start_i = get_cth_index();
    }
}

void cbrtpi_part(int max_size) {
    int start_i = get_cbrtpi_index();

    while (start_i >= 0) {
        int len = start_chunk_size;
        if (start_i + len > max_size) len = max_size - start_i;
        int counter = 0;
        while (counter < len) {
            int k = start_i + counter;
            if (k == 0) {
                m2[k] = cbrtpi(m2[k]);
            } else {
                m2[k] = cbrtpi(m2[k] + m2_copy[k - 1]);
            }
            ++counter;
        }
        start_i = get_cbrtpi_index();
    }
}

void merge_part(int max_size) {
    int start_i = get_merge_index();

    while (start_i >= 0) {
        int len = start_chunk_size;
        if (start_i + len > max_size) len = max_size - start_i;
        int counter = 0;
        while (counter < len) {
            int k = start_i + counter;
            m2[k] = (double) m1[k] * m2[k];
            ++counter;
        }
        start_i = get_merge_index();
    }
}

struct timing time_pr;
void *main_function(void *args) {
    struct main_args thread_args = *((struct main_args *) args);

    int id = thread_args.id;
    unsigned int tmp1 = thread_args.index;
    unsigned int tmp2 = thread_args.index;

    unsigned int local_tmp1 = tmp1;
    unsigned int local_tmp2 = tmp2;
    time_pr.start = get_wtime();
    // GENERATE

    if (id == 0) {
        for (int j = 0; j < N; j++) {
            double value = 1 + rand_r(&local_tmp1) % (A - 1);
            m1[j] = value;
//            printf("%2.f\n", value);
        }
        for (int j = 0; j < N / 2; j++) {
            double value = A + rand_r(&local_tmp2) % (A * 10 - A);
            m2[j] = value;
            m2_copy[j] = value;
//            printf("%.2f\n", value);
        }
    }

    pthread_barrier_wait(&barrier); // join потоков
    time_pr.gen = get_wtime();
    // MAP
    cth_part(N);
    cbrtpi_part(N / 2);

    pthread_barrier_wait(&barrier); // join потоков
    time_pr.map = get_wtime();
    
    // MERGE
    merge_part(N / 2);

    pthread_barrier_wait(&barrier); // join потоков
    time_pr.merge = get_wtime();

    // SORT
    if (id == 0) {
        stupid_sort();
        time_pr.sort = get_wtime();
        double result = 0;
        int j = 0;
        while (j < N / 2 && m2[j] == 0) {
            j++;
        }
        double min = m2[j];
        for (int k = 0; k < N / 2; k++) {
            if (((long) (m2[k] / min) % 2) == 0) {
                result += sin(m2[k]);
            }
        }

//        /*
        pthread_mutex_lock(&print_mutex);
        time_pr.reduce = get_wtime();
        printf("X: %f\n", result);
        pthread_mutex_unlock(&print_mutex);
//         */
    }

    time_pr.end = get_wtime();
    pthread_exit(NULL);
}

void *percent_counter(void *args) {
    int value;
    for (;;) {
        pthread_mutex_lock(&percent_mutex);
        value = *percent;
        printf("Current percent: %d\n", value);
        pthread_mutex_unlock(&percent_mutex);

        if (value >= 100) break;
        sleep(1);
    }
    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    struct timeval T1, T2;
    long delta_ms;

    if (argc < 4) {
        printf("Need to add size of array and number of threads as input arguments\n");
        return -1;
    }

    N = atoi(argv[1]);
    int M = atoi(argv[2]);
    if (M < 1) {
        printf("threads_counter must be positive");
        return -1;
    }

    start_chunk_size = atoi(argv[3]);
    if (start_chunk_size < 1 || start_chunk_size >= N) {
        printf("incorrect start_chuck_size parameter");
        return -2;
    }

    if (argc >= 5) {
        FOR_I = atoi(argv[4]);
    } else FOR_I = 100;

    THREAD_NUM = M;
    if (N < THREAD_NUM) THREAD_NUM = N;

    m1 = (double *) malloc(N * sizeof(double));
    m2 = (double *) malloc(N / 2 * sizeof(double));
    m2_copy = (double *) malloc(N / 2 * sizeof(double));

    percent = (int *) malloc(sizeof(int));
    *percent = 0;

    pthread_mutex_init(&print_mutex, NULL);
    pthread_mutex_init(&chunk_mutex, NULL);
    pthread_barrier_init(&barrier, NULL, THREAD_NUM); //инициализация барьера
    pthread_t percent_tid;
    pthread_t thread[THREAD_NUM];
    struct main_args thread_args[THREAD_NUM];
    pthread_create(&percent_tid, NULL, percent_counter, percent);
    
    gettimeofday(&T1, NULL);

    for (int l = 0; l < FOR_I; l++) {
        cth_index = 0;
        cbrtpi_index = 0;
        merge_index = 0;

        for (int i = 0; i < THREAD_NUM; i++) { // 1, 2 ... THREAD_NUM-1
            thread_args[i].index = l;
            thread_args[i].id = i;

//            printf("create thread #%d\n", i);
            pthread_create(&thread[i], NULL, main_function, &thread_args[i]);
        }
//        printf("main create all\n");


        for (int i = 0; i < THREAD_NUM; i++) {
            pthread_join(thread[i], NULL);
        }
        
        pthread_mutex_lock(&percent_mutex);
        *percent = (100 * (l + 1)) / FOR_I;
        pthread_mutex_unlock(&percent_mutex);
    }

    //Выход из потока:
    pthread_barrier_destroy(&barrier);
 
    pthread_join(percent_tid, NULL);
    gettimeofday(&T2, NULL);
    delta_ms = 1000 * (T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
    print_timing(&time_pr);
    printf("%d;%ld\n", N, delta_ms);
}
