#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <chrono>

struct Point {
    std::vector<float> values;
    int clusterId = -1;
};

struct Centroid {
    std::vector<float> values;
};

float euclideanDistance(const Point& a, const Centroid& b) {
    float sum = 0.0f;
    for (size_t i = 0; i < a.values.size(); ++i)
        sum += (a.values[i] - b.values[i]) * (a.values[i] - b.values[i]);
    return std::sqrt(sum);
}

void updateCentroids(const std::vector<Point>& points, std::vector<Centroid>& centroids, int k, int dim) {
    std::vector<int> counts(k, 0);
    for (auto& centroid : centroids)
        centroid.values.assign(dim, 0.0f);

    for (const auto& p : points) {
        if (p.clusterId >= 0) {
            for (int d = 0; d < dim; ++d)
                centroids[p.clusterId].values[d] += p.values[d];
            counts[p.clusterId]++;
        }
    }

    for (int i = 0; i < k; ++i)
        if (counts[i] > 0)
            for (int d = 0; d < dim; ++d)
                centroids[i].values[d] /= counts[i];
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Uso: " << argv[0] << " <k>\n";
        return 1;
    }

    int k = std::atoi(argv[1]);
    if (k <= 0) {
        std::cerr << "Error: k debe ser un entero positivo.\n";
        return 1;
    }

    std::ifstream file("salida", std::ios::binary);
    if (!file) {
        std::cerr << "No se pudo abrir el archivo de entrada.\n";
        return 1;
    }

    uint32_t nRows, nCols;
    file.read(reinterpret_cast<char*>(&nRows), sizeof(uint32_t));
    file.read(reinterpret_cast<char*>(&nCols), sizeof(uint32_t));

    std::vector<Point> points(nRows);
    for (auto& p : points) {
        p.values.resize(nCols);
        file.read(reinterpret_cast<char*>(p.values.data()), sizeof(float) * nCols);
    }

    file.close();

    const int maxIterations = 2000;
    const float maxDeltaRatio = 0.05f;

    std::vector<Centroid> centroids(k);
    std::srand(std::time(nullptr));
    for (int i = 0; i < k; ++i)
        centroids[i].values = points[i * nRows / k].values;

    for (int i = 0; i < nRows; ++i)
        points[i].clusterId = i % k;

    auto start = std::chrono::high_resolution_clock::now();

    int iteration = 0;
    int changed;
    do {
        changed = 0;
        updateCentroids(points, centroids, k, nCols);

        for (auto& p : points) {
            float minDist = std::numeric_limits<float>::max();
            int bestCluster = -1;
            for (int i = 0; i < k; ++i) {
                float dist = euclideanDistance(p, centroids[i]);
                if (dist < minDist) {
                    minDist = dist;
                    bestCluster = i;
                }
            }
            if (bestCluster != p.clusterId) {
                p.clusterId = bestCluster;
                changed++;
            }
        }
        iteration++;
    } while (changed > nRows * maxDeltaRatio && iteration < maxIterations);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Tiempo de ejecución del algoritmo K-means: " << duration.count() << " segundos\n";
    std::cout << "Terminado en " << iteration << " iteraciones.\n";

    for (int i = 0; i < k; ++i)
        std::cout << "Centroide " << i << ": (" << centroids[i].values[0] << ", " << centroids[i].values[1] << ")\n";

 // ... [código previo intacto]
 std::cout << "\n--- Estadísticas por grupo ---\n";

 for (int c = 0; c < k; ++c) {
     std::vector<float> sum(nCols, 0.0f);
     std::vector<float> sqSum(nCols, 0.0f);
     std::vector<float> minVals(nCols, std::numeric_limits<float>::max());
     std::vector<float> maxVals(nCols, std::numeric_limits<float>::lowest());
     int count = 0;

     for (const auto& p : points) {
         if (p.clusterId == c) {
             for (int d = 0; d < nCols; ++d) {
                 float v = p.values[d];
                 sum[d] += v;
                 sqSum[d] += v * v;
                 minVals[d] = std::min(minVals[d], v);
                 maxVals[d] = std::max(maxVals[d], v);
             }
             count++;
         }
     }

     std::cout << "Grupo " << c << " (" << count << " puntos):\n";
     for (int d = 0; d < nCols; ++d) {
         float media = sum[d] / count;
         float varianza = (sqSum[d] / count) - (media * media);
         std::cout << "  Columna " << d << ": "
                   << "Media=" << media
                   << ", Min=" << minVals[d]
                   << ", Max=" << maxVals[d]
                   << ", Varianza=" << varianza << "\n";
     }
 }

 std::cout << "\n--- Distancia cuadrática media por grupo ---\n";
 for (int c = 0; c < k; ++c) {
     float sumaDistancias = 0.0f;
     int count = 0;
     for (const auto& p : points) {
         if (p.clusterId == c) {
             float dist2 = 0.0f;
             for (int d = 0; d < nCols; ++d) {
                 float diff = p.values[d] - centroids[c].values[d];
                 dist2 += diff * diff;
             }
             sumaDistancias += dist2;
             count++;
         }
     }
     float dcm = (count > 0) ? sumaDistancias / count : 0.0f;
     std::cout << "Grupo " << c << ": DCM = " << dcm << "\n";
 }

 return 0;
}
