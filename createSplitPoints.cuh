#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iterator>
#include <string>


class createSplitPoints{

public:
    createSplitPoints(const point_vector_t& dataset, int NofPartitions, int sampleSize)
        : splits(new point_vector_t){
        
        this->dataset = dataset;
        this->sampleSize = sampleSize;
        this->NofPartitions = NofPartitions;

        this->takeSample();
    
    }

    const point_vector_t& GetSplits() const{
        return *splits;
    }

    int GetNofPartitions() const{
        return NofPartitions;
    }

    void quicksort(point_vector_t& xCoord, int first, int last){

        int i, j, pivot; 
        point temp;

        if( first < last ){
            pivot = first;
            i = first;
            j = last;
      
            while( i < j ) {

                while(xCoord[i].x <= xCoord[pivot].x && i<last)
                    i++;
                while(xCoord[j].x > xCoord[pivot].x)
                    j--;
                if( i < j ){

                  temp = xCoord[i];
                  xCoord[i] = xCoord[j];
                  xCoord[j] = temp;
               
                }
            }
      
            temp = xCoord[pivot];
            xCoord[pivot] = xCoord[j];
            xCoord[j] = temp;
            quicksort(xCoord, first, j-1);
            quicksort(xCoord, j+1, last);
      
        }

    }

    void mapPartitionsToThreads(point_vector_t& dataset, int* indexesOfSplits, point_vector_t& splitsPoints, int NofPartitions, int NofPoints){

        int i, j;

        i = 0;
        for(j=0; j < NofPartitions-1; j++){
        
            while (dataset[i].x != splitsPoints[j].x){
                i++;
            }   
            indexesOfSplits[j] = i;
        }
        //the index of the last split point will be the one with the biggest x-coordinate
        indexesOfSplits[j] = NofPoints-1;

    }

protected:
    int sampleSize;
    int NofPartitions;
    point_vector_t dataset;
    point_vector_t sampledDataset;

    /** \brief Template method for loading data files. It can add the points in an internal or external memory vector
     *
     * \param dataset the dataset from which we take sample and we want to partition
     * \param SplitDataset SplitPointVector& the internal or external memory vector to add split points into
     * \param sampleSize the size of the sample we want for finding the split points
     */
     template<class SplitPointVector>
     void sampleTheDataset(const point_vector_t& dataset, point_vector_t& sampledDataset, SplitPointVector& SplitDataset, int sampleSize, int NofPartitions)
     {
        
        int randomi;
        sampledDataset.reserve(sampleSize);
        for (int i=0; i < sampleSize; i++){
            
            randomi = rand()%sampleSize + 1;
            sampledDataset[i] = dataset[randomi];
            
        }

        findSplitPoints(sampledDataset, SplitDataset, sampleSize, NofPartitions);

     }
   

private:

    std::unique_ptr<point_vector_t> splits;

    void takeSample(){

        sampleTheDataset(dataset, sampledDataset, *splits, sampleSize, NofPartitions);
         
    }

    template<class SplitPointVector>
    void findSplitPoints(point_vector_t& sampledDataset, SplitPointVector& SplitDataset, int sampleSize, int NofPartitions){

        int step;

        quicksort(sampledDataset, 0, sampleSize-1);

        step = (int) sampleSize / NofPartitions;

        SplitDataset.reserve(NofPartitions);

        for(int j=1; j < NofPartitions; j++){

            SplitDataset.push_back(sampledDataset[step*j]);
        }
        
    }

};