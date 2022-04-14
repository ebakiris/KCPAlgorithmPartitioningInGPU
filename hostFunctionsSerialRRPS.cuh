void check(int error){
    if(error == 0){
        printf("successful memory allocation!\n");
    }
    if(error == 2){
        printf("The API call failed because it was unable to allocate enough memory to perform the requested operation.\n");
    }
    if(error == 1){
        printf("This indicates that one or more of the parameters passed to the API call is not within an acceptable range of values.\n");
    }
    if(error == 21){
        printf("This indicates that the direction of the memcpy passed to the API call is not one of the types specified by cudaMemcpyKind.\n");
    }
}

double EuclideanHost(double x1, double y1, double x2, double y2)
{
	double x = x1 - x2; //calculating number to square in next step
	double y = y1 - y2;
	double dist;

	dist = pow(x, 2) + pow(y, 2);       //calculating Euclidean distance
	dist = sqrt(dist);                  

	return dist;
}

bool isFullHost(double* localMaxKHeap, int NofCP){

    int i;
    bool value = true;
    
    
    for(i = 0; i < NofCP; i++){    
        if(localMaxKHeap[i] <= DBL_MIN){
            value = false;
            break;
        }
    }      
    
    return value;

}

double maxDelta(double* localMaxKHeap, int NofCP){

    double max = localMaxKHeap[0];

    for(int i=1; i < NofCP; i++){

        if(localMaxKHeap[i] > max){
            max = localMaxKHeap[i];
        }
    }

    return max;

}

int insertDistanceHost(double dist, double* localMaxKHeap, bool isfull, int NofCP){

    double max;
    int i;
    double epsilon = DBL_MIN; 
    
    if(isfull == false){

        i = 0;
        while(!(localMaxKHeap[i] <= epsilon)){
            i++;
        }
        localMaxKHeap[i] = dist;
        return i;
        
    }
    else{
        
        
        int index = 0;
        max = localMaxKHeap[0];

        for(i=1; i < NofCP; i++){

            if(localMaxKHeap[i] > max){
                max = localMaxKHeap[i];
                index = i;
            }
        }

        if(dist < max & dist > epsilon){

            localMaxKHeap[index] = dist;
            return index;
        }


    }
    return -1;
}

void quicksort(double* xCoord, int first, int last){

    int i, j, pivot; 
    double temp;
    

    if( first < last ){
        pivot = first;
        i = first;
        j = last;
  
        while( i < j ) {

            while(xCoord[i] <= xCoord[pivot] && i<last)
                i++;
            while(xCoord[j] > xCoord[pivot])
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

void quicksortVectors(point_vector_t& xCoord, int first, int last){

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
        quicksortVectors(xCoord, first, j-1);
        quicksortVectors(xCoord, j+1, last);
  
    }
}

void quicksortPairs(closestpairs* xCoord, int first, int last){

    int i, j, pivot; 
    closestpairs temp;

    if( first < last ){
        pivot = first;
        i = first;
        j = last;
  
        while( i < j ) {

            while(xCoord[i].dist <= xCoord[pivot].dist && i<last)
                i++;
            while(xCoord[j].dist > xCoord[pivot].dist)
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
        quicksortPairs(xCoord, first, j-1);
        quicksortPairs(xCoord, j+1, last);
  
    }

}


void deleteSamePairs(closestpairs* globalMaxKHeap, int sizeOfArray){

    int i;
    int j;


    for(i = 0; i < sizeOfArray; i++){

        for(j = i+1; j < sizeOfArray; j++){

            if(globalMaxKHeap[i].p.id == globalMaxKHeap[j].p.id & globalMaxKHeap[i].q.id == globalMaxKHeap[j].q.id){

                globalMaxKHeap[j].dist = DBL_MIN;
            }

        }
    }
    
}

void deleteSameDists(double* globalMaxKHeap, int sizeOfArray){

    int i;
    int j;


    for(i = 0; i < sizeOfArray; i++){

        for(j = i+1; j < sizeOfArray; j++){

            if(globalMaxKHeap[i] == globalMaxKHeap[j]){

                globalMaxKHeap[j] = DBL_MIN;
            }

        }
    }
    
}
