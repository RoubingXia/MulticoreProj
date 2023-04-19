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

void write_helper1(char file_name[], double alpha_betas[9801][2], double pointers[9801][14], int best_alpha_beta_idx) {
    // write an array to a file, len1: length of the Smoothed value for month,
    // len2: total length of each double array, len3: total length of the data array
    FILE *fptr;
    fptr = fopen(file_name,"w");
    if(fptr == NULL)
    {
        printf("Error!");
        exit(1);
    }
    // pointers is an array of double pointers, alpha_betas stores pointers which points to double arrays where alpha and beta are stored. iterate the data

        //double data[] = pointers[i];    // get data
        //double info[] = alpha_betas[i]; // get corresponding alpha and beta
        double alpha = alpha_betas[best_alpha_beta_idx][0];
        double beta = alpha_betas[best_alpha_beta_idx][1];
        fprintf(fptr,"Alpha : %f Beta : %f\n",alpha, beta);
        for (int j = 0; j < 14; ++j) {
            int month = j + 1;
            double line =  pointers[best_alpha_beta_idx][j];
            if (j >= 12) {
                // Output the forecasted values
                fprintf(fptr,"\tPredicted sales for month %d: %.2lf\n",month, line);
            }
            else {
                // Output the smoothed values for the first 11 months
                fprintf(fptr,"\tSmoothed value for month %d: %.2lf\n",month , line);
            }

        }



    fclose(fptr);
}


void write_helper(char file_name[], double alpha_betas[9801][2], double pointers[9801][14]) {
    // write an array to a file, len1: length of the Smoothed value for month,
    // len2: total length of each double array, len3: total length of the data array
    FILE *fptr;
    fptr = fopen(file_name,"w");
    if(fptr == NULL)
    {
        printf("Error!");
        exit(1);
    }
    // pointers is an array of double pointers, alpha_betas stores pointers which points to double arrays where alpha and beta are stored. iterate the data
    for (int i = 0; i < 9801; ++i) {
        //double data[] = pointers[i];    // get data
        //double info[] = alpha_betas[i]; // get corresponding alpha and beta
        double alpha = alpha_betas[i][0];
        double beta = alpha_betas[i][1];
        fprintf(fptr,"Alpha : %f Beta : %f\n",alpha, beta);
        for (int j = 0; j < 14; ++j) {
            int month = j + 1;
            double line =  pointers[i][j];
            if (j >= 11) {
                // Output the forecasted values
                fprintf(fptr,"\tPredicted sales for month %d: %.2lf\n",month, line);
            }
            else {
                // Output the smoothed values for the first 11 months
                fprintf(fptr,"\tSmoothed value for month %d: %.2lf\n",month , line);
            }

        }
    }


    fclose(fptr);
}

double getMean(double* source_data, double* smooth_data){
    double mean = 0.0;
    for (int i = 0; i < 10; ++i) {
        mean += abs(source_data[i] - smooth_data[i]);
    }
    return mean / 10;
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
    int best_alpha_beta_idx = 0;
    double MAE = 100000.0;
    double alpha = 0.01;
    double beta = 0.01;
    double base = 0.01;
    double t_start = 0.0, t_taken;

    int range = 100;
    double data_set[(range - 1) * (range - 1)][11 + 3]; // output data
    double alpha_betas[(range - 1) * (range - 1)][2]; // alpha beta
    int idx = 0;
    for (int i = 1; i < range; ++i) {
        for (int j = 1; j < range; ++j) {
            alpha_betas[idx][0] = base * i;
            alpha_betas[idx++][1] = base * j;
        }
    }
    int len = (range - 1) * (range - 1); // length of alpha_betas
    t_start = omp_get_wtime();

    if (threads_count == 0) {
        int idx = 0;// used to fill pointers and alpha_betas
        for (int i = 0 ; i < len; ++i) {
            alpha = alpha_betas[i][0];
            beta = alpha_betas[i][1];
            double_exponential_smoothing(data, data_length, alpha, beta, data_set[i]);
            // get current mean
            double mean = getMean(data, data_set[i]);
            if (mean < MAE) {
                MAE = mean;
                best_alpha_beta_idx = i;
            }
        }
    }
    else {
        int step = len / threads_count;
        #pragma omp parallel num_threads(threads_count)
        {
            int tid = omp_get_thread_num();
            int start = tid * step;
            int end = (tid == threads_count - 1) ? len : (tid + 1) * step;
            for (int i = start ; i < end; ++i) {
                alpha = alpha_betas[i][0];
                beta = alpha_betas[i][1];
                double_exponential_smoothing(data, data_length, alpha, beta, data_set[i]);
                double mean = getMean(data, data_set[i]);
                #pragma omp critical
                {
                    if (mean < MAE) {
                        MAE = mean;
                        best_alpha_beta_idx = i;
                    }
                }

            }
        }
    }


    t_taken =  omp_get_wtime() - t_start;
    printf("Time taken : %f \n", t_taken);
    /*
double result[14];
    double_exponential_smoothing(data, data_length, 0.01, 0.03, result);
    // Output the smoothed values for the first 11 months
    for (int i = 0; i < data_length; i++) {
        printf("Smoothed value for month %d: %.2lf\n", i + 1, result[i]);
    }

    // Output the forecasted values
    for (int i = 0; i < 3; i++) {
        printf("Predicted sales for month %d: %.2lf\n", data_length + 1 + i, result[data_length + i]);
    }
     */
    write_helper("output", alpha_betas, data_set);
    write_helper1("optimal_output", alpha_betas, data_set, best_alpha_beta_idx);
    printf("best alpha : %.2lf beta : %.2lf\n", alpha_betas[best_alpha_beta_idx][0], alpha_betas[best_alpha_beta_idx][1]);
    return 0;
}
