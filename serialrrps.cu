#include <iostream>
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <type_traits>
#include <memory>
#include <array>
#include <chrono>
#include <cfloat>
#include <time.h>
#include "KCPProblem.cuh"
#include "hostFunctionsSerialRRPS.cuh"
#include "rrpsAlgorithm.cuh"
using namespace std::chrono; 


int main(int argc, char* argv[]){

    //we need the number of closest pairs (k) to be specified from command line
    if ( argc < 6 ) {
        
        std::cout << "Argument error. Please enter:\n";
        std::cout << "Argument 1: The number of the closest pairs:\n";
        std::cout << "Argument 2: The file of the first dataset P:\n";
        std::cout << "Argument 3: The file of the second dataset Q:\n";
        std::cout << "Argument 4: The number of the total points of dataset P\n";
        std::cout << "Argument 5: The number of the total points of dataset Q\n";
        
        
        return 1;
    }

    try{

        int NofCP = atoi(argv[1]);
        int NofPointsP = atoi(argv[4]);
        int NofPointsQ = atoi(argv[5]);
        
        point_vector_t datasetP, datasetQ;
        point* h_datasetP = (point*)malloc((NofPointsP+1)*sizeof(point));
        point* h_datasetQ = (point*)malloc((NofPointsQ+1)*sizeof(point));
        double* globalMaxKHeap = (double*)malloc(NofCP*sizeof(double));
    
        
        // Default values
        std::string fileDatasetP("");
        std::string fileDatasetQ("");

        fileDatasetP = argv[2];
        fileDatasetQ = argv[3];

        std::unique_ptr<KCPProblem> query;
        
        //allocate the problem object depending
        //read the files and store the datasets 
        query.reset(new KCPProblem(fileDatasetP, fileDatasetQ, NofCP));
    
        //start time
        auto start = high_resolution_clock::now();

        datasetP = query->GetDatasetP();
        datasetQ = query->GetDatasetQ();
        
        // sort the two datasets
        quicksortVectors(datasetP, 0, NofPointsP - 1);
        quicksortVectors(datasetQ, 0, NofPointsQ - 1);
        
        
        for(int i=0; i < NofPointsP; i++){


            h_datasetP[i].id = datasetP[i].id;
            h_datasetP[i].x = datasetP[i].x;
            h_datasetP[i].y = datasetP[i].y;
        }

        for(int i=0; i < NofPointsQ; i++){

            h_datasetQ[i].id = datasetQ[i].id;
            h_datasetQ[i].x = datasetQ[i].x;
            h_datasetQ[i].y = datasetQ[i].y;

        }

        h_datasetP[NofPointsP].id = NofPointsP;
        h_datasetP[NofPointsP].x = 0.0;
        h_datasetP[NofPointsP].y = 0.0;

        h_datasetQ[NofPointsQ].id = NofPointsQ;
        h_datasetQ[NofPointsQ].x = 0.0;
        h_datasetQ[NofPointsQ].y = 0.0;
        
        
        ReverseRunPlaneSweepSerial(h_datasetP, h_datasetQ, NofPointsP, NofPointsQ, NofCP, globalMaxKHeap);
        
        quicksort(globalMaxKHeap, 0, NofCP-1);

        for(int i=0; i < NofCP ; i++){
            printf("%.10lf\n", globalMaxKHeap[i]);
        }
        
        auto stop = high_resolution_clock::now();

        std::cout << "Total Execution Time (Serial RRPS)= " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "[μs]" << std::endl;
      
    }
    catch (std::exception & ex)
    {
        //report any exception
        std::cout << "Exception: " << ex.what() << std::endl;
        return 1;
    }


}
