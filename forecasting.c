#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

int read_csv(const char *filename, double *data) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file: %s\n", filename);
        exit(1);
    }

    // Read and discard the header row
    fscanf(file, "%*[^\n]\n");

    int i = 0;
    while (fscanf(file, "%*d,%lf", &data[i]) != EOF) {
        printf("Read value: %lf\n", data[i]);
        i++;
    }
    fclose(file);
    return i;
}

void double_exponential_smoothing(const double *data, int n, double alpha, double beta, double *result) {
    double lt[n];
    double bt[n];

    // Initialize the level and trend values
    lt[0] = data[0];
    bt[0] = data[1]-data[0];

    // Calculate the level and trend values for the time series
    for (int i = 1; i < n; i++) {
        lt[i] = alpha * data[i] + (1 - alpha) * (lt[i - 1] + bt[i - 1]);
        bt[i] = beta * (lt[i] - lt[i - 1]) + (1 - beta) * bt[i - 1];
    }

    // Store the smoothed values in the result array
    for (int i = 0; i < n; i++) {
        result[i] = lt[i] + bt[i];
    }

    // Predict the values for the next 3 months
    int max_month = 3;
    for (int h = 1; h <= max_month; h++) {
        result[n - 1 + h] = lt[n - 1] + h * bt[n - 1];
    }
}



int main(int argc, char *argv[]) {

    if(argc != 2 ){
        printf("usage: ./forecasting t N t\n");
        printf("t: the number of threads \n");
        exit(1);
    }
    int threads_count = atoi(argv[1]);
    const char *filename = "output.csv";
    double data[11];
    int data_length = read_csv(filename, data);

    double alpha = 0.01;
    double beta = 0.01;
    double t_start = 0.0, t_taken;
    double result[11 + 3];
    int range = 99;

    t_start = omp_get_wtime();

    if (threads_count == 0) {
        for (int f_a = 1; f_a < 100; ++f_a) {
            for (int f_b = 1; f_b < 100; ++f_b) {
                alpha *= f_a;
                beta *= f_b;
                double_exponential_smoothing(data, data_length, alpha, beta, result);
            }
        }
    }
    else {
        int step = range / threads_count;
        #pragma omp parallel num_threads(threads_count)
        {

            int tid = omp_get_thread_num();
            int start = tid * step;
            int end = (tid == threads_count - 1) ? range - 1 : (tid + 1) * step;
            for (int f_a = start; f_a <= end; ++f_a) {
                for (int f_b = start; f_b <= end; ++f_b) {
                    alpha *= f_a;
                    beta *= f_b;
                    double_exponential_smoothing(data, data_length, alpha, beta, result);
                }
            }
        }
    }


    t_taken =  omp_get_wtime() - t_start;
    printf("Time taken : %f \n", t_taken);
    // Output the smoothed values for the first 11 months
    for (int i = 0; i < data_length; i++) {
        printf("Smoothed value for month %d: %.2lf\n", i + 1, result[i]);
    }

    // Output the forecasted values
    for (int i = 0; i < 3; i++) {
        printf("Predicted sales for month %d: %.2lf\n", data_length + 1 + i, result[data_length + i]);
    }

    return 0;
}