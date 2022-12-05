#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <fwBase.h>
#include <fwSignal.h>

//глупая сортировка
int stupid_sort(float arr[], int size) {
   int i;
   for(i = 0; i < size - 1; i++) {
      if(arr[i] > arr[i + 1]){
         int tmp = arr[i];
         arr[i] = arr[i + 1];
         arr[i + 1] = tmp;
         i = -1;
      }
   }
   return 0;
}

//гиперболический котангенс корня числа
float cth(int x){
    return (exp(2*sqrt(x))+1)/(exp(2*sqrt(x))-1);
}

//Кубический корень после умножения на число Пи
float cbrtpi(int x){
    return cbrt(x * M_PI);
}

float min_but_not_zero(float arr[], int size) {
   int i = 0;
   float min;
   do {
      min = arr[i];
      i++;
   } while(min == 0 && i < size);

   //all elements are equal to zero
   if(min == 0){
      return 0;
   }

   for(i = 0; i < size; i++){
      if(arr[i] != 0 && arr[i] < min)
         min = arr[i];
   }
   return min;
}

float reduce(float arr[], int size) {
   float sum = 0;
   int i;
   float min = min_but_not_zero(arr, size);

   if(min == 0){
      return -1;
   }
   
   for(i = 0; i < size; i++){
      if(((int)(arr[i]/min)) % 2 == 0)
         sum += sin(arr[i]);
   }
   return sum;
}

int main(int argc, char* argv[])
{
    int i, N, M;
    struct timeval T1, T2;
    long delta_ms;
    unsigned int seeds = 320;

    N = atoi(argv[1]); /* N равен первому параметру командной строки */
    M = atoi(argv[2]); /* M - второй аргумент, количество потоков*/
    gettimeofday(&T1, NULL); /* запомнить текущее время T1 */
    fwSetNumThreads(M);
    //N = 1160;
    float M1[N];
    float M2[N/2];
    float M2_copy[N/2];
    
    for (i=0; i<100; i++) /* 100 экспериментов */
    {
        srand(i); /* инициализировать начальное значение ГСЧ */
        /* Заполнить массив исходных данных размером N */

        for(int m1=0; m1<N; m1++){
            M1[m1] = rand_r(&seeds)%320+1;       
        }

        for(int m2=0; m2<N/2; m2++){
            M2[m2] = rand_r(&seeds)%2881+320;
            M2_copy[m2] = M2[m2];
            //printf("INIT VALUE: m1: %f, m1_copy: %f\n", M2[m2], M2_copy[m2]);
        }

        /* Решить поставленную задачу, заполнить массив с результатами */
        //Array M1
        // for(int m1=0; m1<N; m1++){
        //     M1[m1] = cth(M1[m1]);
        // }

        //cth wit fw functions
        fwsSqrt_32f_I(M1, N);
        fwsTanh_32f_A24(M1, M1, N);
        fwsDivCRev_32f_I(1, M1, N);

        //ArrayM2
        // for(int m2=1; m2<N/2; m2++){
        //     M2[m2] = M2[m2] + M2_copy[m2-1];
        //     printf("plus result: %f\n", M2[m2]);
        // }

        fwsAdd_32f(M2+1, M2_copy, M2+1, N/2);

        // for(int m2=1; m2<N/2; m2++){
        //     M2[m2] = M2[m2] + M2_copy[m2-1];
        // }


        // for(int m2=0; m2<N/2; m2++){
        //     M2[m2] = cbrtpi(M2[m2]);
        // }

        fwsMulC_32f(M2, M_PI, M2, N/2);
        fwsCubrt_32f(M2, M2, N/2);

        // //Array M1[m]*M2[2]
        // for(int m2=0; m2<N/2; m2++){
        //     M2_copy2[m2] = M1[m2]*M2[m2];
        // }
        
        fwsMul_32f(M1, M2, M2, N/2);

        /* Отсортировать массив с результатами указанным методом */\
        stupid_sort(M2, N/2);
    }
    gettimeofday(&T2, NULL); /* запомнить текущее время T2 */
    delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;
    printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms); /* T2 - T1 */
    printf("\nVerification=%f.\n", reduce(M2, N/2));
    
    return 0;
}
