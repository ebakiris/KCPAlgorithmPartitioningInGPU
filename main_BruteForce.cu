#include <cstdlib>
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
#include "createSplitPoints.cuh"
#include "deviceFunctions.cuh"
#include "hostFunctions.cuh"
#include "kernelPhase2.cuh"
#include "kernelSample.cuh"
#include "BruteForceAlgorithm.cuh"
using namespace std::chrono; 

#define THREADS_PER_BLOCK 256

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
        
        point_vector_t splitsP, splitsQ;
        point_vector_t datasetP, datasetQ;
        
        point* h_datasetP = (point*)malloc(NofPointsP*sizeof(point));
        point* h_datasetQ = (point*)malloc(NofPointsQ*sizeof(point));
       
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
        auto startBrute = high_resolution_clock::now(); 
        
        datasetP = query->GetDatasetP();
        datasetQ = query->GetDatasetQ();
        
        if(NofPointsP == NofPointsQ){
            for(int i=0; i < NofPointsP; i++){
                h_datasetP[i] = datasetP[i];            
                h_datasetQ[i] = datasetQ[i];
            }
        }
        else{
            for(int i=0; i < NofPointsP; i++){
                h_datasetP[i] = datasetP[i];            
            }
            for(int i=0; i < NofPointsQ; i++){
                h_datasetQ[i] = datasetQ[i];
            }
        }

    
        BruteForceAlgorithm(h_datasetP, h_datasetQ, NofCP, NofPointsP, NofPointsQ);
        auto stopBrute = high_resolution_clock::now();

        std::cout << "Time difference (Brute Force) = " << std::chrono::duration_cast<std::chrono::milliseconds>(stopBrute - startBrute).count() << "[ms]" << std::endl;
        
        
        
    }
    catch (std::exception & ex)
    {
        //report any exception
        std::cout << "Exception: " << ex.what() << std::endl;
        return 1;
    }


}
