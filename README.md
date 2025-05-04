
# Práctica 1 – Algoritmo K-means (CAP)

Este proyecto contiene tres implementaciones del algoritmo K-means:

- **Versión Secuencial** (`algoritmo-normal.cpp`)
- **Versión paralelizada con OpenMP** (`algoritmo-OpenMP.cpp`)
- **Versión paralelizada con OpenMPI** (`algoritmo-OpenMPI.cpp`)

Además, se incluye un generador de datos de entrada en binario (`generateRndSeeds.cpp`).

---

## 📦 Contenido

- `generateRndSeeds.cpp`: Generador de conjuntos de datos en formato binario (`salida`)
- `algoritmo-normal.cpp`: Implementación secuencial
- `algoritmo-OpenMP.cpp`: Implementación paralela con OpenMP
- `algoritmo-OpenMPI.cpp`: Implementación paralela con OpenMPI
- `salida`: Archivo binario con los datos generados

---

## ⚙️ Requisitos

- **C++11 o superior**
- **OpenMP** (para `algoritmo-OpenMP.cpp`)
- **OpenMPI** (para `algoritmo-OpenMPI.cpp`)

---

## 🛠️ Compilación

### 1. Generador de datos

```bash
g++ -std=c++11 -o generateRnd generateRndSeeds.cpp
./generateRnd
```

Esto generará un archivo binario `salida` con los datos.

---

### 2. Versión Secuencial

```bash
g++ -std=c++11 -o kmeans_normal algoritmo-normal.cpp
./kmeans_normal 4
```

---

### 3. Versión con OpenMP

```bash
g++ -fopenmp -std=c++11 -o kmeans_omp algoritmo-OpenMP.cpp
./kmeans_omp 4
```

---

### 4. Versión con OpenMPI

```bash
mpic++ -std=c++11 -o kmeans_mpi algoritmo-OpenMPI.cpp
mpirun -np 4 ./kmeans_mpi
```

---

## 📊 Salida esperada

Cada implementación mostrará:

- Tiempo de ejecución
- Iteraciones realizadas
- Coordenadas de los centroides
- Estadísticas por grupo: media, mínimo, máximo y varianza por columna
- DCM (distancia cuadrática media) por grupo

---

## 📁 Entregables (según PDF de la práctica)

1. Memoria técnica (formato PDF)
2. Archivo comprimido con:
   - `algoritmo-normal.cpp`
   - `algoritmo-OpenMP.cpp`
   - `algoritmo-OpenMPI.cpp`
   - `generateRndSeeds.cpp`
   - `README.md`
   - (opcional) `salida` de ejemplo

---

## 🧠 Notas

- El valor de `k` debe coincidir con el número de procesos o hilos deseado (normalmente 4).
- El archivo de entrada debe llamarse exactamente `salida` y estar en formato binario como especifica el enunciado.
- Se puede volver a generar con `generateRnd`.

