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
        std::cout << "Argument 4: The number of the total partitions\n";
        std::cout << "Argument 5: The number of the total points of dataset P\n";
        std::cout << "Argument 6: The number of the total points of dataset Q\n";
        
        
        return 1;
    }

    try{

        int NofCP = atoi(argv[1]);
        int NofPartitions = atoi(argv[4]);
        int NofPointsP = atoi(argv[5]);
        int NofPointsQ = atoi(argv[6]);
        //NofSamplePartitions is actually the sample from which we will calculate the bound in phase1(kernelSample)
        int NofSamplePartitions = (int) NofPartitions*0.1;
        //SampleSize is for creating the split points(partitions)
        int SampleSizeP = (int) NofPointsP*0.3;
        int SampleSizeQ = (int) NofPointsQ*0.3;
        point_vector_t splitsP, splitsQ;
        point_vector_t datasetP, datasetQ;
        pair_vector_t neccArray ;
        point* h_datasetP = (point*)malloc(NofPointsP*sizeof(point));
        point* h_datasetQ = (point*)malloc(NofPointsQ*sizeof(point));
        //the maximum size of the necessary pairs of partitions to check 
        pair* h_neccArray = (pair*)malloc(NofPartitions*NofPartitions*sizeof(pair));
        ovpairs* neccPairs = (ovpairs*)malloc(NofSamplePartitions*sizeof(ovpairs));

        int *indexesOfSplitsP = new int[NofPartitions];
        int *indexesOfSplitsQ = new int[NofPartitions];
        int *indexesOfSplitsPOfSample = new int[NofSamplePartitions];
        int *indexesOfSplitsQOfSample = new int[NofSamplePartitions];
        
        // Default values
        std::string fileDatasetP("");
        std::string fileDatasetQ("");

        fileDatasetP = argv[2];
        fileDatasetQ = argv[3];

        std::unique_ptr<KCPProblem> query;
        
        //allocate the problem object depending
        //read the files and store the datasets 
        query.reset(new KCPProblem(fileDatasetP, fileDatasetQ, NofCP));
        
        std::unique_ptr<createSplitPoints> PointsP;
        std::unique_ptr<createSplitPoints> PointsQ;

        auto startTotal = high_resolution_clock::now();

        //find the splits points for each dataset P and Q
        PointsP.reset(new createSplitPoints(query->GetDatasetP(), NofPartitions, SampleSizeP));
        PointsQ.reset(new createSplitPoints(query->GetDatasetQ(), NofPartitions, SampleSizeQ));
        
        splitsP = PointsP->GetSplits();
        splitsQ = PointsQ->GetSplits();
        
        // sort the two datasets
        datasetP = query->GetDatasetP();
        PointsP->quicksort(datasetP, 0, NofPointsP-1);
        
        datasetQ = query->GetDatasetQ();
        PointsQ->quicksort(datasetQ, 0, NofPointsQ-1);
        
        //map the indexes of split points of partitions to each thread
        PointsP->mapPartitionsToThreads(datasetP, indexesOfSplitsP, splitsP, NofPartitions, NofPointsP);
        PointsQ->mapPartitionsToThreads(datasetQ, indexesOfSplitsQ, splitsQ, NofPartitions, NofPointsQ);
        
        // Invoke kernel for sample
        int threadsPerBlock = NofSamplePartitions;
        int blocksPerGrid = 1;
        int sizeofheap = threadsPerBlock*blocksPerGrid; 

        double* globalMaxKHeap = (double*)malloc(sizeofheap*NofCP*sizeof(double));
        

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

        for(int j=0; j < NofSamplePartitions; j++){
            indexesOfSplitsPOfSample[j] = indexesOfSplitsP[j];
            indexesOfSplitsQOfSample[j] = indexesOfSplitsQ[j];
            neccPairs[j].priority = -1;
        }

        // Allocate vectors in device memory
        double* d_globalMaxKHeap;
        point* d_datasetP;
        point* d_datasetQ;
        pair* d_neccArray;
        closestpairs* phase2_d_globalMaxKHeap;
        ovpairs* d_neccPairs;
        
        
        //start time of bound calculation
        auto start = high_resolution_clock::now(); 

        findOverlapPartitions(datasetP, datasetQ, indexesOfSplitsP, indexesOfSplitsQ, neccPairs, NofSamplePartitions, NofPartitions);

        cudaMalloc(&d_globalMaxKHeap, sizeofheap*NofCP*sizeof(double));
        cudaMalloc(&d_datasetP, NofPointsP*sizeof(point));
        cudaMalloc(&d_datasetQ, NofPointsQ*sizeof(point));
        cudaMalloc(&d_neccPairs, NofSamplePartitions*sizeof(ovpairs));
        
        // Copy vectors from host memory to device memory
        cudaMemcpy(d_datasetP, h_datasetP, NofPointsP*sizeof(point), cudaMemcpyHostToDevice);
        cudaMemcpy(d_datasetQ, h_datasetQ, NofPointsQ*sizeof(point), cudaMemcpyHostToDevice);
        cudaMemcpy(d_neccPairs, neccPairs, NofSamplePartitions*sizeof(ovpairs), cudaMemcpyHostToDevice);

        // Invoke kernel
        //start time of sample phase
        auto startPhaseSample = high_resolution_clock::now();
        kernelSample<<<blocksPerGrid, threadsPerBlock>>>(d_datasetP, d_datasetQ, NofCP, d_globalMaxKHeap, d_neccPairs);
        cudaDeviceSynchronize();
        auto stopPhaseSample = high_resolution_clock::now();
        
        std::cout << "CLEAN Time difference(Sample phase) = " << std::chrono::duration_cast<std::chrono::milliseconds>(stopPhaseSample - startPhaseSample).count() << "[ms]" << std::endl;
        
        // Copy result from device memory to host memory
        cudaMemcpy(globalMaxKHeap, d_globalMaxKHeap, sizeofheap*NofCP*sizeof(double), cudaMemcpyDeviceToHost);

        auto stop = high_resolution_clock::now();
        //end time of bound calculation

        std::cout << "Time difference(Sample for bound) = " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "[ms]" << std::endl;

        quicksort(globalMaxKHeap, 0, sizeofheap*NofCP-1);

        int count = 0;
        int g;
        for(g=0; g < sizeofheap*NofCP; g++){
            if(globalMaxKHeap[g] != DBL_MIN){
                count++;
            }
            if(count == NofCP){
                break;
            }
        }

        printf("The bound found by sample is: %.10f\n", globalMaxKHeap[g]);

        double bound = globalMaxKHeap[g];
        
        cudaFree(d_globalMaxKHeap);
        cudaFree(d_neccPairs);

        neccArray = findNeccPartitions(datasetP, datasetQ, indexesOfSplitsP, indexesOfSplitsQ, bound, neccArray, NofPartitions);
       
        printf("The total pair of partitions for calculations are: %d\n", neccArray.size());

        for(int i = 0; i < neccArray.size(); i++){
            h_neccArray[i] = neccArray[i];
        }

        //determine the number of blocks and number of threads depending on the 
        //neccessary pairs for calculations, that we extract based on bound
        int notUsefulThreadsId;
        int numberOfPairs = neccArray.size();
        if(numberOfPairs > THREADS_PER_BLOCK){
            threadsPerBlock = THREADS_PER_BLOCK;
            if(numberOfPairs % THREADS_PER_BLOCK != 0){
                blocksPerGrid = (numberOfPairs / THREADS_PER_BLOCK) + 1; 
                notUsefulThreadsId = numberOfPairs;
            }
            else{
                blocksPerGrid = (numberOfPairs / THREADS_PER_BLOCK);
                notUsefulThreadsId = -1;
            }
        }
        else{
            threadsPerBlock = numberOfPairs;
            blocksPerGrid = 1;
            notUsefulThreadsId = -1;
        }

        closestpairs* phase2_globalMaxKHeap = (closestpairs*)malloc(numberOfPairs*NofCP*sizeof(closestpairs));

        //start time of phase2
        auto startPhase2 = high_resolution_clock::now(); 

        cudaMalloc(&d_neccArray, numberOfPairs*sizeof(pair));
        cudaMalloc(&phase2_d_globalMaxKHeap, numberOfPairs*NofCP*sizeof(closestpairs));
        // Copy vectors from host memory to device memory
        cudaMemcpy(d_neccArray, h_neccArray, numberOfPairs*sizeof(pair), cudaMemcpyHostToDevice);
        
        //start time of phase2
        auto startPhase2Clean = high_resolution_clock::now(); 
        kernelPhase2<<<blocksPerGrid, threadsPerBlock>>>(d_datasetP, d_datasetQ, d_neccArray, NofCP, phase2_d_globalMaxKHeap, notUsefulThreadsId, bound, numberOfPairs);
        cudaDeviceSynchronize();
        auto stopPhase2Clean = high_resolution_clock::now();

        std::cout << "CLEAN Time difference(Phase2 final) = " << std::chrono::duration_cast<std::chrono::milliseconds>(stopPhase2Clean - startPhase2Clean).count() << "[ms]" << std::endl;
        
        cudaMemcpy(phase2_globalMaxKHeap, phase2_d_globalMaxKHeap, numberOfPairs*NofCP*sizeof(closestpairs), cudaMemcpyDeviceToHost);
        
        auto stopPhase2 = high_resolution_clock::now();

        std::cout << "Time difference(Phase2 final) = " << std::chrono::duration_cast<std::chrono::milliseconds>(stopPhase2 - startPhase2).count() << "[ms]" << std::endl;
        //end time of phase2
        
        auto startSort = high_resolution_clock::now();

        quicksortPairs(phase2_globalMaxKHeap, 0, numberOfPairs*NofCP-1);
        //insertionSort(phase2_globalMaxKHeap, numberOfPairs*NofCP);
        //bubbleSort(phase2_globalMaxKHeap, numberOfPairs*NofCP);
        //selectionSort(phase2_globalMaxKHeap, numberOfPairs*NofCP);
        //mergeSort(phase2_globalMaxKHeap, 0, numberOfPairs*NofCP-1);

        auto stopSort = high_resolution_clock::now();

        std::cout << "Sort time = " << std::chrono::duration_cast<std::chrono::milliseconds>(stopSort - startSort).count() << "[ms]" << std::endl;

        //deleteSamePairs(phase2_globalMaxKHeap, numberOfPairs*NofCP);

        int i = 0;
        int counter = 0;
        while(counter < NofCP){
            if(phase2_globalMaxKHeap[i].dist != DBL_MIN){

                if(i != 0){
                    if(phase2_globalMaxKHeap[i].p.id != phase2_globalMaxKHeap[i-1].p.id || phase2_globalMaxKHeap[i].q.id != phase2_globalMaxKHeap[i-1].q.id){
                        counter++;
                        printf("The pair is: p (%f, %f) and q (%f, %f) with Distance: %.10lf\n", phase2_globalMaxKHeap[i].p.x, phase2_globalMaxKHeap[i].p.y, phase2_globalMaxKHeap[i].q.x, phase2_globalMaxKHeap[i].q.y, phase2_globalMaxKHeap[i].dist);
                    }    
                }
                else{
                    counter++;
                    printf("The pair is: p (%f, %f) and q (%f, %f) with Distance: %.10lf\n", phase2_globalMaxKHeap[i].p.x, phase2_globalMaxKHeap[i].p.y, phase2_globalMaxKHeap[i].q.x, phase2_globalMaxKHeap[i].q.y, phase2_globalMaxKHeap[i].dist);
                }
            
            }
            i++;            
        }

        // Free device memory
        cudaFree(d_datasetQ);
        cudaFree(d_datasetP);
        cudaFree(d_neccArray);
        cudaFree(phase2_d_globalMaxKHeap);
        
        auto stopTotal = high_resolution_clock::now();

        std::cout << "Total Execution Time = " << std::chrono::duration_cast<std::chrono::milliseconds>(stopTotal - startTotal).count() << "[ms]" << std::endl;
        
        
    }
    catch (std::exception & ex)
    {
        //report any exception
        std::cout << "Exception: " << ex.what() << std::endl;
        return 1;
    }


}
