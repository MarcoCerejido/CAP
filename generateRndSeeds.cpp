#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstdio>

#define PI 3.14159265359f

// Estructura genérica para un punto multidimensional
typedef std::vector<float> Point;

Point getRandomPoint(const Point& centroide, float maxRadius, float minRadius = 0.0f) {
    size_t dim = centroide.size();
    Point p(dim);
    
    float r = minRadius + (maxRadius - minRadius) * ((float)rand() / RAND_MAX);
    float alpha = 2.0f * PI * ((float)rand() / RAND_MAX);

    // Para 2D usamos coordenadas polares, para >2D perturbamos cada dimensión
    if (dim == 2) {
        p[0] = centroide[0] + r * cos(alpha);
        p[1] = centroide[1] + r * sin(alpha);
    } else {
        for (size_t i = 0; i < dim; ++i)
            p[i] = centroide[i] + ((float)rand() / RAND_MAX - 0.5f) * 2.0f * maxRadius;
    }

    return p;
}

int main() {
    srand((unsigned int)time(nullptr));

    int nClusters, nPointsPerCluster, nCol;
    std::cout << "Número de columnas (dimensiones): ";
    std::cin >> nCol;
    std::cout << "Número de clusters: ";
    std::cin >> nClusters;
    std::cout << "Número de puntos por cluster: ";
    std::cin >> nPointsPerCluster;

    std::vector<Point> data;

    for (int i = 0; i < nClusters; ++i) {
        Point centroide(nCol);
        for (int d = 0; d < nCol; ++d)
            centroide[d] = -10.0f + 20.0f * ((float)rand() / RAND_MAX);

        for (int j = 0; j < nPointsPerCluster; ++j)
            data.push_back(getRandomPoint(centroide, 1.0f));
    }

    FILE* resultsFile = fopen("salida", "wb");
    if (!resultsFile) {
        std::cerr << "No se pudo abrir el archivo 'salida' para escribir." << std::endl;
        return 1;
    }

    int nFilas = nClusters * nPointsPerCluster;
    fwrite(&nFilas, sizeof(int), 1, resultsFile);
    fwrite(&nCol, sizeof(int), 1, resultsFile);

    for (const auto& p : data)
        fwrite(p.data(), sizeof(float), nCol, resultsFile);

    fclose(resultsFile);

    for (const auto& p : data) {
        for (float val : p)
            std::cout << val << "\t";
        std::cout << "\n";
    }

    std::cout << "\nGenerados " << nFilas << " puntos en " << nClusters << " clústeres de " << nCol << " dimensiones.\n";
    return 0;
}
