#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <omp.h>

using namespace std;

struct Punto {
    vector<float> valores;
    int grupo;

    Punto(int dimensiones) : valores(dimensiones), grupo(-1) {}
};

float distancia(const Punto& p, const vector<float>& centroide) {
    float suma = 0.0f;
    for (size_t i = 0; i < p.valores.size(); ++i)
        suma += pow(p.valores[i] - centroide[i], 2);
    return sqrt(suma);
}

vector<float> calcularCentroide(const vector<Punto>& puntos, int grupo, int dimensiones) {
    vector<float> centroide(dimensiones, 0.0f);
    int total = 0;

    for (const auto& p : puntos) {
        if (p.grupo == grupo) {
            for (int i = 0; i < dimensiones; ++i)
                centroide[i] += p.valores[i];
            total++;
        }
    }

    if (total > 0) {
        for (int i = 0; i < dimensiones; ++i)
            centroide[i] /= total;
    }

    return centroide;
}

void leerArchivo(const string& filename, vector<Punto>& puntos, int& dimensiones) {
    ifstream archivo(filename, ios::binary);
    if (!archivo) {
        cerr << "Error al abrir el archivo " << filename << endl;
        exit(1);
    }

    int filas;
    archivo.read(reinterpret_cast<char*>(&filas), sizeof(int));
    archivo.read(reinterpret_cast<char*>(&dimensiones), sizeof(int));
    puntos.resize(filas, Punto(dimensiones));

    for (int i = 0; i < filas; ++i) {
        archivo.read(reinterpret_cast<char*>(puntos[i].valores.data()), dimensiones * sizeof(float));
    }

    archivo.close();
}

void kmeans(vector<Punto>& puntos, int k, int dimensiones) {
    const int max_iter = 2000;
    const float max_desplazados_ratio = 0.05f;
    int total_puntos = puntos.size();
    vector<vector<float>> centroides(k, vector<float>(dimensiones));

    srand(time(nullptr));
    for (int i = 0; i < k; ++i) {
        int idx = rand() % total_puntos;
        centroides[i] = puntos[idx].valores;
    }

    bool continuar = true;
    int iter = 0;
    double start = omp_get_wtime();

    while (continuar && iter < max_iter) {
        int desplazados = 0;

        #pragma omp parallel for reduction(+:desplazados)
        for (int i = 0; i < total_puntos; ++i) {
            float min_dist = distancia(puntos[i], centroides[0]);
            int mejor_grupo = 0;

            for (int j = 1; j < k; ++j) {
                float dist = distancia(puntos[i], centroides[j]);
                if (dist < min_dist) {
                    min_dist = dist;
                    mejor_grupo = j;
                }
            }

            if (puntos[i].grupo != mejor_grupo) {
                desplazados++;
                puntos[i].grupo = mejor_grupo;
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < k; ++j) {
            centroides[j] = calcularCentroide(puntos, j, dimensiones);
        }

        float desplazado_ratio = float(desplazados) / total_puntos;
        continuar = desplazado_ratio >= max_desplazados_ratio;
        iter++;
    }

    double end = omp_get_wtime();
    double tiempo = end - start;

    cout << "\nTiempo de ejecución del algoritmo K-means: " << tiempo << " segundos" << endl;
    cout << "Terminado en " << iter << " iteraciones.\n";

    for (int i = 0; i < k; ++i) {
        cout << "Centroide " << i << ": (";
        for (int d = 0; d < dimensiones; ++d) {
            cout << centroides[i][d];
            if (d < dimensiones - 1) cout << ", ";
        }
        cout << ")\n";
    }

    cout << "\n--- Estadísticas por grupo ---\n";

    for (int grupo = 0; grupo < k; ++grupo) {
        vector<float> suma(dimensiones, 0.0f), min_vals(dimensiones, INFINITY), max_vals(dimensiones, -INFINITY), varianza(dimensiones, 0.0f);
        int count = 0;

        for (const auto& p : puntos) {
            if (p.grupo != grupo) continue;
            count++;
            for (int d = 0; d < dimensiones; ++d) {
                float val = p.valores[d];
                suma[d] += val;
                min_vals[d] = min(min_vals[d], val);
                max_vals[d] = max(max_vals[d], val);
            }
        }

        vector<float> media(dimensiones);
        for (int d = 0; d < dimensiones; ++d)
            media[d] = (count > 0) ? suma[d] / count : 0.0f;

        for (const auto& p : puntos) {
            if (p.grupo != grupo) continue;
            for (int d = 0; d < dimensiones; ++d) {
                float diff = p.valores[d] - media[d];
                varianza[d] += diff * diff;
            }
        }

        for (int d = 0; d < dimensiones; ++d)
            varianza[d] = (count > 0) ? varianza[d] / count : 0.0f;

        cout << "Grupo " << grupo << " (" << count << " puntos):\n";
        for (int d = 0; d < dimensiones; ++d) {
            cout << "  Columna " << d << ": Media=" << media[d]
                 << ", Min=" << min_vals[d]
                 << ", Max=" << max_vals[d]
                 << ", Varianza=" << varianza[d] << endl;
        }
    }

    cout << "\n--- Distancia cuadrática media por grupo ---\n";
    for (int grupo = 0; grupo < k; ++grupo) {
        float sumaDist2 = 0.0f;
        int count = 0;
        for (const auto& p : puntos) {
            if (p.grupo == grupo) {
                float dist2 = 0.0f;
                for (int d = 0; d < dimensiones; ++d) {
                    float diff = p.valores[d] - centroides[grupo][d];
                    dist2 += diff * diff;
                }
                sumaDist2 += dist2;
                count++;
            }
        }
        float dcm = (count > 0) ? sumaDist2 / count : 0.0f;
        cout << "Grupo " << grupo << ": DCM = " << dcm << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Uso: " << argv[0] << " <k>\n";
        return 1;
    }

    int k = atoi(argv[1]);
    if (k <= 0) {
        cerr << "Error: k debe ser un entero positivo.\n";
        return 1;
    }

    vector<Punto> puntos;
    int dimensiones;
    leerArchivo("salida", puntos, dimensiones);  // Archivo fijo
    kmeans(puntos, k, dimensiones);

    return 0;
}
