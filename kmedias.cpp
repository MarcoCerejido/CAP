//este es el base de chat gpt


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ITER 2000
#define THRESHOLD 0.05 // 5% cambio

typedef struct {
    float *data;
    int rows, cols;
} Dataset;

typedef struct {
    float *centroids;
    int k;
} KMeans;

void read_binary_file(const char *filename, Dataset *dataset) {
    FILE *file = fopen(filename, "rb");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    
    fread(&dataset->rows, sizeof(unsigned int), 1, file);
    fread(&dataset->cols, sizeof(unsigned int), 1, file);
    
    dataset->data = (float *)malloc(dataset->rows * dataset->cols * sizeof(float));
    fread(dataset->data, sizeof(float), dataset->rows * dataset->cols, file);
    
    fclose(file);
}

void initialize_kmeans(KMeans *kmeans, Dataset *dataset, int k) {
    kmeans->k = k;
    kmeans->centroids = (float *)malloc(k * dataset->cols * sizeof(float));
    
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < dataset->cols; j++) {
            kmeans->centroids[i * dataset->cols + j] = dataset->data[i * dataset->cols + j];
        }
    }
}

void update_centroids(KMeans *kmeans, Dataset *dataset, int *assignments) {
    float *new_centroids = (float *)calloc(kmeans->k * dataset->cols, sizeof(float));
    int *counts = (int *)calloc(kmeans->k, sizeof(int));
    
    for (int i = 0; i < dataset->rows; i++) {
        int cluster = assignments[i];
        for (int j = 0; j < dataset->cols; j++) {
            new_centroids[cluster * dataset->cols + j] += dataset->data[i * dataset->cols + j];
        }
        counts[cluster]++;
    }
    
    for (int i = 0; i < kmeans->k; i++) {
        for (int j = 0; j < dataset->cols; j++) {
            if (counts[i] > 0) {
                kmeans->centroids[i * dataset->cols + j] = new_centroids[i * dataset->cols + j] / counts[i];
            }
        }
    }
    
    free(new_centroids);
    free(counts);
}

void compute_statistics(Dataset *dataset) {
    int cols = dataset->cols;
    float min[cols], max[cols], mean[cols], variance[cols];
    
    for (int j = 0; j < cols; j++) {
        min[j] = max[j] = dataset->data[j];
        mean[j] = 0.0;
    }
    
    for (int i = 0; i < dataset->rows; i++) {
        for (int j = 0; j < cols; j++) {
            float val = dataset->data[i * cols + j];
            if (val < min[j]) min[j] = val;
            if (val > max[j]) max[j] = val;
            mean[j] += val;
        }
    }
    
    for (int j = 0; j < cols; j++) {
        mean[j] /= dataset->rows;
    }
    
    for (int i = 0; i < dataset->rows; i++) {
        for (int j = 0; j < cols; j++) {
            float val = dataset->data[i * cols + j] - mean[j];
            variance[j] += val * val;
        }
    }
    
    for (int j = 0; j < cols; j++) {
        variance[j] /= dataset->rows;
    }
    
    printf("Statistics:\n");
    for (int j = 0; j < cols; j++) {
        printf("Col %d -> Min: %.2f, Max: %.2f, Mean: %.2f, Var: %.2f\n", j, min[j], max[j], mean[j], variance[j]);
    }
}

int main(int argc, char **argv) {
    Dataset dataset;
    read_binary_file("data.bin", &dataset);
    
    KMeans kmeans;
    initialize_kmeans(&kmeans, &dataset, 3);
    
    int *assignments = (int *)malloc(dataset.rows * sizeof(int));
    for (int i = 0; i < MAX_ITER; i++) {
        update_centroids(&kmeans, &dataset, assignments);
    }
    
    compute_statistics(&dataset);
    
    free(dataset.data);
    free(kmeans.centroids);
    free(assignments);
    
    return 0;
}
