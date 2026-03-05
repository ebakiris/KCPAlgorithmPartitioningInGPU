# Parallel Processing of Spatial Queries (K-Closest Pairs)

A high-performance implementation of spatial query processing using parallel computing architectures. This project focuses on solving the **K-Closest Pairs Query (K-CPQ)** by leveraging the massive parallelism of NVIDIA GPUs via the CUDA platform.

## 📖 Overview

Spatial queries are essential in fields such as urban planning, traffic routing, and environmental modeling. This project implements a parallel algorithm designed to efficiently find the $K$ pairs of points with the smallest Euclidean distances between two large datasets.

By partitioning data and utilizing a specialized **"Reverse Run Plane-Sweep" (RRPS)** technique, the algorithm significantly reduces the computational overhead compared to traditional brute-force methods.

## 🚀 Key Features

* **GPU Acceleration:** Optimized for NVIDIA architectures to handle large-scale spatial datasets.


* **Advanced Plane-Sweep:** Implements an improved "Reverse Run Plane-Sweep" with a "sliding window" to minimize distance calculations.


* **Data Partitioning:** Efficiently divides datasets into spatial partitions to maximize thread utilization and load balancing.


* **Scalability:** Tested with synthetic and real-world datasets (from SpatialHadoop) containing up to 1 million points.



## 🛠️ Tech Stack

* **Language:** C / C++
* **Parallel Computing:** NVIDIA CUDA 


* **Environment:** Linux (tested on NVIDIA Quadro P400) 



## 🔧 Installation & Setup

### Prerequisites

* A system running **Linux**.


* An **NVIDIA GPU** with CUDA support.


* **CUDA Toolkit** installed.



### Compilation

You can compile the project using `nvcc`:

```bash
nvcc -o main main.cu 

```

*(Note: Replace `main.cu` with your primary source file name if different.)*

## 💻 Usage

The program requires several parameters to define the query and dataset constraints:

```bash
./main <K> <File_P> <File_Q> <Partitions> <Size_P> <Size_Q>

```

**Example:**
To find the 10 closest pairs between two datasets of 500,000 points each using 2,850 partitions:

```bash
./main 10 P_gaussian_500K.txt Q_gaussian_500K.txt 2850 500000 500000

```



## 📊 Performance

The project includes a detailed experimental evaluation comparing serial and parallel executions. For massive datasets, the parallel RRPS algorithm provides a substantial speedup over brute-force (BF) approaches.

## 🎓 Background

This project was developed as a Diploma Thesis at the **University of Thessaly**, Department of Electrical and Computer Engineering.

* **Author:** Emmanouil Bakiris 


* **Supervisor:** Michael Vassilakopoulos 


* **Year:** 2021


The full research can be found here: https://www.e-ce.uth.gr/wp-content/uploads/formidable/59/Bakiris_emmanouil.pdf

## 📜 License & Copyright

### Code License
The source code is licensed under the **MIT License**.


### Academic Rights & Attribution
In accordance with the internal regulations of the University of Thessaly:

* **Academic Credit:** Any use, reproduction, or distribution of this work for academic or research purposes must provide full attribution to the author (Emmanouil Bakiris) and the supervisor (Michael Vassilakopoulos), citing the original thesis title: "Parallel Processing of Spatial Queries".

* **Institutional Ownership:** Intellectual property rights are jointly held by the author and the Department of Electrical and Computer Engineering at the University of Thessaly.

* **Commercial Use:** For any commercial exploitation or use of the algorithms/codebase, written permission must be obtained from both the author and the supervising professor.



