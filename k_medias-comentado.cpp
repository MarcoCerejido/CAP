// //esta es la versión """"""""""""mejorada"""""""""""" por mi



// #include <stdio.h>
// #include <stdlib.h>
// #define MAX_ITER 2000

// typedef struct {
//     float *data;
//     int rows, cols;
// } Dataset;

// typedef struct {
//     float *centroids;
//     int k;
// } KMeans;

// void leer_fichero_binario(const char *nombre, Dataset *dataset) {
//     FILE *file = fopen(nombre, "rb");
//     if (!file) {
//         perror("No se encuentra el archivo");
//         exit(EXIT_FAILURE);
//     }

//     // Leer número de filas y columnas
//     fread(&dataset->rows, sizeof(unsigned int), 1, file);
//     fread(&dataset->cols, sizeof(unsigned int), 1, file);
 

//     // Asignar memoria para los datos
//     dataset->data = (float *)malloc(dataset->rows * dataset->cols * sizeof(float));

//     // Leer los datos de la matriz
//     fread(dataset->data, sizeof(float), dataset->rows * dataset->cols, file);
//     fclose(file);
// }



// void ini_kmeans(KMeans *kmeans, Dataset *dataset, int k) {
//     kmeans->k = k;
//     kmeans->centroids = (float *)malloc(k * dataset->cols * sizeof(float));

//     // Inicializar los centroides con los primeros k puntos del dataset
//     for (int i = 0; i < k; i++) {
//         for (int j = 0; j < dataset->cols; j++) {
//             kmeans->centroids[i * dataset->cols + j] = dataset->data[i * dataset->cols + j];
//         }
//     }
// }

// void update_centroids(KMeans *kmeans, Dataset *dataset, int *assignments) {
//     float *new_centroids = (float *)calloc(kmeans->k * dataset->cols, sizeof(float));
//     int *counts = (int *)calloc(kmeans->k, sizeof(int));

//     // Asignar cada punto a su centroide más cercano
//     for (int i = 0; i < dataset->rows; i++) {
//         float min_distance = __FLT_MAX__;
//         int cluster = -1;
//         for (int j = 0; j < kmeans->k; j++) {
//             float distance = 0;
//             for (int d = 0; d < dataset->cols; d++) {
//                 float diff = dataset->data[i * dataset->cols + d] - kmeans->centroids[j * dataset->cols + d];
//                 distance += diff * diff;
//             }
//             if (distance < min_distance) {
//                 min_distance = distance;
//                 cluster = j;
//             }
//         }
//         assignments[i] = cluster;
//         for (int j = 0; j < dataset->cols; j++) {
//             new_centroids[cluster * dataset->cols + j] += dataset->data[i * dataset->cols + j];
//         }
//         counts[cluster]++;
//     }

//     // Recalcular los centroides
//     for (int i = 0; i < kmeans->k; i++) {
//         for (int j = 0; j < dataset->cols; j++) {
//             if (counts[i] > 0) {
//                 kmeans->centroids[i * dataset->cols + j] = new_centroids[i * dataset->cols + j] / counts[i];
//             }
//         }
//     }

//     free(new_centroids);
//     free(counts);
// }

// void compute_statistics(Dataset *dataset) {
//     int cols = dataset->cols;
//     float min[cols], max[cols], mean[cols], variance[cols];

//     // Inicializar los valores de las estadísticas
//     for (int j = 0; j < cols; j++) {
//         min[j] = max[j] = dataset->data[j];
//         mean[j] = 0.0;
//         variance[j] = 0.0;
//     }

//     // Calcular mínimo, máximo y media
//     for (int i = 0; i < dataset->rows; i++) {
//         for (int j = 0; j < cols; j++) {
//             float val = dataset->data[i * cols + j];
//             if (val < min[j]) min[j] = val;
//             if (val > max[j]) max[j] = val;
//             mean[j] += val;
//         }
//     }

//     for (int j = 0; j < cols; j++) {
//         mean[j] /= dataset->rows;
//     }

//     // Calcular varianza
//     for (int i = 0; i < dataset->rows; i++) {
//         for (int j = 0; j < cols; j++) {
//             float val = dataset->data[i * cols + j] - mean[j];
//             variance[j] += val * val;
//         }
//     }

//     for (int j = 0; j < cols; j++) {
//         variance[j] /= dataset->rows;
//     }

//     // Imprimir las estadísticas
//     printf("Statistics:\n");
//     for (int j = 0; j < cols; j++) {
//         printf("Col %d -> Min: %.2f, Max: %.2f, Mean: %.2f, Var: %.2f\n", j, min[j], max[j], mean[j], variance[j]);
//     }
// }

// int main() {
//     Dataset dataset;
//     leer_fichero_binario("data.bin", &dataset);

//     KMeans kmeans;
//     ini_kmeans(&kmeans, &dataset, 3);  // Usamos k = 3 como ejemplo

//     int *assignments = (int *)malloc(dataset.rows * sizeof(int));
    
//     // Ejecutar el algoritmo k-medias
//     for (int i = 0; i < MAX_ITER; i++) {
//         update_centroids(&kmeans, &dataset, assignments);
//     }

//     // Calcular estadísticas
//     compute_statistics(&dataset);

//     free(dataset.data);
//     free(kmeans.centroids);
//     free(assignments);

//     return 0;
// }
