
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <limits>

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

    if (rank == 0) std::cout << "K-means con OpenMPI usando " << size << " procesos.\n";

    int nCols = 0;
    std::vector<Point> localPoints;
    std::vector<Centroid> centroids(size);

    // Nodo 0 lee los datos y los distribuye
    if (rank == 0) {
        std::ifstream file("salida", std::ios::binary);
        if (!file) {
            std::cerr << "Error al abrir el archivo.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int nRows;
        file.read(reinterpret_cast<char*>(&nRows), sizeof(int));
        file.read(reinterpret_cast<char*>(&nCols), sizeof(int));

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
        int localRows;
        MPI_Recv(&localRows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Bcast(&nCols, 1, MPI_INT, 0, MPI_COMM_WORLD);

        localPoints.resize(localRows);
        std::vector<float> buffer(localRows * nCols);
        MPI_Recv(buffer.data(), localRows * nCols, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < localRows; ++i)
            localPoints[i].values = std::vector<float>(buffer.begin() + i * nCols, buffer.begin() + (i + 1) * nCols);
    }

    MPI_Bcast(&nCols, 1, MPI_INT, 0, MPI_COMM_WORLD);

    for (auto& p : localPoints)
        p.clusterId = rank;

    const int maxIter = 2000;
    const float minRatio = 0.05f;
    int iter = 0;

    bool converged = false;
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

        int totalMoved = 0, totalLocal = localPoints.size();
        MPI_Allreduce(&moved, &totalMoved, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        int totalPoints = 0;
        MPI_Allreduce(&totalLocal, &totalPoints, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        float ratio = static_cast<float>(totalMoved) / totalPoints;
        converged = (ratio < minRatio);
        iter++;
    }

    if (rank == 0) std::cout << "\nConvergencia tras " << iter << " iteraciones.\n";

    // Cálculo de DCM local
    Centroid myCentroid;
    computeCentroid(localPoints, myCentroid, nCols, rank);

    float sumDist2 = 0.0f;
    int count = 0;
    for (const auto& p : localPoints) {
        if (p.clusterId == rank) {
            sumDist2 += distanceSquared(p, myCentroid);
            count++;
        }
    }

    float dcm = (count > 0) ? sumDist2 / count : 0.0f;
    std::cout << "Proceso " << rank << " - DCM = " << dcm << ", puntos = " << count << "\n";

    MPI_Finalize();
    return 0;
}
