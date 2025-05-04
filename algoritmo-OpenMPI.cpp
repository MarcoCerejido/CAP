#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <iomanip>

struct Point {
    std::vector<float> values;
    int clusterId = -1;
};

using Centroid = std::vector<float>;

float distanceSquared(const Point& p, const Centroid& c) {
    float dist = 0.0f;
    for (size_t i = 0; i < p.values.size(); ++i) {
        float diff = p.values[i] - c[i];
        dist += diff * diff;
    }
    return dist;
}

void computeCentroid(const std::vector<Point>& points, Centroid& centroid, int dim, int clusterId) {
    centroid.assign(dim, 0.0f);
    int count = 0;

    for (const auto& p : points) {
        if (p.clusterId == clusterId) {
            for (int d = 0; d < dim; ++d)
                centroid[d] += p.values[d];
            count++;
        }
    }

    if (count > 0) {
        for (int d = 0; d < dim; ++d)
            centroid[d] /= count;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int nCols = 0, localRows = 0;
    std::vector<Point> localPoints;

    if (rank == 0) {
        std::ifstream file("salida", std::ios::binary);
        if (!file) {
            std::cerr << "Error al abrir el archivo." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int nRows;
        file.read(reinterpret_cast<char*>(&nRows), sizeof(int));
        file.read(reinterpret_cast<char*>(&nCols), sizeof(int));
        MPI_Bcast(&nCols, 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::vector<float> allData(nRows * nCols);
        for (int i = 0; i < nRows; ++i)
            file.read(reinterpret_cast<char*>(&allData[i * nCols]), sizeof(float) * nCols);
        file.close();

        int blockSize = nRows / size;
        for (int i = 0; i < size; ++i) {
            int start = i * blockSize;
            int count = (i == size - 1) ? (nRows - start) : blockSize;
            int nSend = count * nCols;
            if (i == 0) {
                localRows = count;
                localPoints.resize(count);
                for (int j = 0; j < count; ++j)
                    localPoints[j].values = std::vector<float>(allData.begin() + (start + j) * nCols,
                                                               allData.begin() + (start + j + 1) * nCols);
            } else {
                MPI_Send(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(allData.data() + start * nCols, nSend, MPI_FLOAT, i, 1, MPI_COMM_WORLD);
            }
        }
    } else {
        MPI_Bcast(&nCols, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Recv(&localRows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        localPoints.resize(localRows);
        std::vector<float> buffer(localRows * nCols);
        MPI_Recv(buffer.data(), localRows * nCols, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < localRows; ++i)
            localPoints[i].values = std::vector<float>(buffer.begin() + i * nCols, buffer.begin() + (i + 1) * nCols);
    }

    for (auto& p : localPoints)
        p.clusterId = rank;

    const int maxIter = 2000;
    const float minRatio = 0.05f;
    int iter = 0;
    bool converged = false;
    double startTime = MPI_Wtime();

    while (!converged && iter < maxIter) {
        Centroid localCentroid;
        computeCentroid(localPoints, localCentroid, nCols, rank);

        std::vector<float> allCentroids(nCols * size);
        MPI_Allgather(localCentroid.data(), nCols, MPI_FLOAT, allCentroids.data(), nCols, MPI_FLOAT, MPI_COMM_WORLD);

        std::vector<Centroid> centroids(size, std::vector<float>(nCols));
        for (int i = 0; i < size; ++i)
            for (int d = 0; d < nCols; ++d)
                centroids[i][d] = allCentroids[i * nCols + d];

        int moved = 0;
        for (auto& p : localPoints) {
            float minDist = std::numeric_limits<float>::max();
            int bestId = -1;
            for (int i = 0; i < size; ++i) {
                float d = distanceSquared(p, centroids[i]);
                if (d < minDist) {
                    minDist = d;
                    bestId = i;
                }
            }
            if (bestId != p.clusterId) {
                p.clusterId = bestId;
                moved++;
            }
        }

        int totalMoved, totalPoints = localPoints.size(), globalPoints;
        MPI_Allreduce(&moved, &totalMoved, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&totalPoints, &globalPoints, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        float ratio = static_cast<float>(totalMoved) / globalPoints;
        converged = (ratio < minRatio);
        iter++;
    }

    double endTime = MPI_Wtime();
    double elapsed = endTime - startTime;

    // Estadísticas por grupo
    std::vector<float> stats(4 * nCols, 0.0f);  // media, min, max, varianza por columna
    std::vector<float> mins(nCols, std::numeric_limits<float>::max());
    std::vector<float> maxs(nCols, std::numeric_limits<float>::lowest());
    std::vector<float> sums(nCols, 0.0f);
    std::vector<float> sqsums(nCols, 0.0f);

    for (const auto& p : localPoints) {
        if (p.clusterId == rank) {
            for (int d = 0; d < nCols; ++d) {
                float val = p.values[d];
                sums[d] += val;
                sqsums[d] += val * val;
                if (val < mins[d]) mins[d] = val;
                if (val > maxs[d]) maxs[d] = val;
            }
        }
    }

    int count = 0;
    for (const auto& p : localPoints)
        if (p.clusterId == rank) count++;

    Centroid myCentroid;
    computeCentroid(localPoints, myCentroid, nCols, rank);

    float dcm = 0.0f;
    for (const auto& p : localPoints)
        if (p.clusterId == rank) dcm += distanceSquared(p, myCentroid);
    dcm = (count > 0) ? dcm / count : 0.0f;

    // Recolectar todo en el proceso 0
    std::vector<float> allCentroids(nCols * size);
    std::vector<float> allSums(nCols * size), allSqSums(nCols * size);
    std::vector<float> allMins(nCols * size), allMaxs(nCols * size);
    std::vector<int> allCounts(size);
    std::vector<float> allDCMs(size);

    MPI_Gather(myCentroid.data(), nCols, MPI_FLOAT, allCentroids.data(), nCols, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(sums.data(), nCols, MPI_FLOAT, allSums.data(), nCols, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(sqsums.data(), nCols, MPI_FLOAT, allSqSums.data(), nCols, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(mins.data(), nCols, MPI_FLOAT, allMins.data(), nCols, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(maxs.data(), nCols, MPI_FLOAT, allMaxs.data(), nCols, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&count, 1, MPI_INT, allCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&dcm, 1, MPI_FLOAT, allDCMs.data(), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Tiempo de ejecución del algoritmo K-means: " << elapsed << " segundos" << std::endl;
        std::cout << "Terminado en " << iter << " iteraciones." << std::endl;
        for (int i = 0; i < size; ++i) {
            std::cout << "Centroide " << i << ": (";
            for (int d = 0; d < nCols; ++d)
                std::cout << allCentroids[i * nCols + d] << (d == nCols - 1 ? "" : ", ");
            std::cout << ")" << std::endl;
        }

        std::cout << "\n--- Estadísticas por grupo ---" << std::endl;
        for (int i = 0; i < size; ++i) {
            std::cout << "Grupo " << i << " (" << allCounts[i] << " puntos):" << std::endl;
            for (int d = 0; d < nCols; ++d) {
                float mean = allSums[i * nCols + d] / allCounts[i];
                float var = (allSqSums[i * nCols + d] / allCounts[i]) - (mean * mean);
                std::cout << "  Columna " << d << ": Media=" << mean
                          << ", Min=" << allMins[i * nCols + d]
                          << ", Max=" << allMaxs[i * nCols + d]
                          << ", Varianza=" << var << std::endl;
            }
        }

        std::cout << "\n--- Distancia cuadrática media por grupo ---" << std::endl;
        for (int i = 0; i < size; ++i)
            std::cout << "Grupo " << i << ": DCM = " << allDCMs[i] << std::endl;
    }

    MPI_Finalize();
    return 0;
}
