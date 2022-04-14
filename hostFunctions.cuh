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

bool isFullHost(void* localMaxKHeap, int NofCP){

    int i;
    bool value = true;
    
    closestpairs* local = (closestpairs*) localMaxKHeap;
    for(i = 0; i < NofCP; i++){    
        if(local[i].dist <= DBL_MIN){
            value = false;
            break;
        }
    }      
    
    return value;

}

int insertDistanceHost(double dist, void* localMaxKHeap, bool isfull, int NofCP){

    double max;
    int i;
    double epsilon = DBL_MIN; 
    
    if(isfull == false){

        
        closestpairs* local = (closestpairs*) localMaxKHeap;
        i = 0;
        while(!(local[i].dist <= epsilon)){
            i++;
        }
        local[i].dist = dist;
        return i;
        
    }
    else{
        
        
        closestpairs* local = (closestpairs*) localMaxKHeap;
        
        int index = 0;
        max = local[0].dist;

        for(i=1; i < NofCP; i++){

            if(local[i].dist > max){
                max = local[i].dist;
                index = i;
            }
        }

        if(dist < max & dist > epsilon){

            local[index].dist = dist;
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

void insertionSort(closestpairs arr[], int n){  
    
    closestpairs key;
    int i, j;  

    for (i = 1; i < n; i++) 
    {  
        key = arr[i];  
        j = i - 1;  
  
        /* Move elements of arr[0..i-1], that are  
        greater than key, to one position ahead  
        of their current position */
        while (j >= 0 && arr[j].dist > key.dist) 
        {  
            arr[j + 1] = arr[j];  
            j = j - 1;  
        }  
        arr[j + 1] = key;  
    }  
}


void swap(closestpairs *xp, closestpairs *yp)  
{  
    closestpairs temp = *xp;  
    *xp = *yp;  
    *yp = temp;  
}  
  
// A function to implement bubble sort  
void bubbleSort(closestpairs arr[], int n)  
{  
    int i, j;  
    for (i = 0; i < n-1; i++)      
      
    // Last i elements are already in place  
    for (j = 0; j < n-i-1; j++)  
        if (arr[j].dist > arr[j+1].dist)  
            swap(&arr[j], &arr[j+1]);  
}

void selectionSort(closestpairs arr[], int n)  
{  
    int i, j, min_idx;  
  
    // One by one move boundary of unsorted subarray  
    for (i = 0; i < n-1; i++)  
    {  
        // Find the minimum element in unsorted array  
        min_idx = i;  
        for (j = i+1; j < n; j++)  
        if (arr[j].dist < arr[min_idx].dist)  
            min_idx = j;  
  
        // Swap the found minimum element with the first element  
        swap(&arr[min_idx], &arr[i]);  
    }  
}

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(closestpairs arr[], int l, int m, int r)
{
    int n1 = m - l + 1;
    int n2 = r - m;
 
    // Create temp arrays
    closestpairs L[n1], R[n2];
 
    // Copy data to temp arrays L[] and R[]
    for (int i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
 
    // Merge the temp arrays back into arr[l..r]
 
    // Initial index of first subarray
    int i = 0;
 
    // Initial index of second subarray
    int j = 0;
 
    // Initial index of merged subarray
    int k = l;
 
    while (i < n1 && j < n2) {
        if (L[i].dist <= R[j].dist) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
 
    // Copy the remaining elements of
    // L[], if there are any
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }
 
    // Copy the remaining elements of
    // R[], if there are any
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}
 
// l is for left index and r is
// right index of the sub-array
// of arr to be sorted */
void mergeSort(closestpairs arr[],int l,int r){
    if(l>=r){
        return;//returns recursively
    }
    int m = (l+r-1)/2;
    mergeSort(arr,l,m);
    mergeSort(arr,m+1,r);
    merge(arr,l,m,r);
}



pair_vector_t findNeccPartitions(point_vector_t datasetP, point_vector_t datasetQ, int* indexesOfSplitsP, int* indexesOfSplitsQ, double bound, pair_vector_t neccArray, int NofPartitions){

    double dist;
    int i;
    int j; 
    int count = 0;
    bool rule1 = false;
    bool rule2 = false;
    bool rule3 = false;
    bool rule4 = false;
    bool readyTobreak = false;
    

    for(i=0; i < NofPartitions; i++){

        for(j=i; j < NofPartitions; j++){

            if(i == 0 ){
                
                //the two partitions overlap
                if(datasetQ[indexesOfSplitsQ[j]].x <= datasetP[indexesOfSplitsP[i]].x & datasetQ[indexesOfSplitsQ[j]].x >= datasetP[0].x){
                    rule1 = true;
                }
              
            }
            else{
                //the two partitions overlap
                if(datasetQ[indexesOfSplitsQ[j]].x <= datasetP[indexesOfSplitsP[i]].x & datasetQ[indexesOfSplitsQ[j]].x >= datasetP[indexesOfSplitsP[i-1]].x ){
                    rule1 = true;
                }

                
            }

            if(rule1 == false){

                //the partitions overlap
                if(j != 0){
                    if(i != 0){
                        if(datasetQ[indexesOfSplitsQ[j-1]].x <= datasetP[indexesOfSplitsP[i]].x & datasetQ[indexesOfSplitsQ[j-1]].x >= datasetP[indexesOfSplitsP[i-1]].x ){
                            rule2 = true;
                        }
                    }
                    else{
                        if(datasetQ[indexesOfSplitsQ[j-1]].x <= datasetP[indexesOfSplitsP[i]].x & datasetQ[indexesOfSplitsQ[j-1]].x >= datasetP[0].x ){
                            rule2 = true;
                        }
                    }
                    
                }
                else{
                    if(i == 0){
                        
                        if(datasetQ[0].x <= datasetP[indexesOfSplitsP[i]].x & datasetQ[0].x >= datasetP[0].x ){
                            rule2 = true;
                        }
                    }
                    else{
                        
                        if(datasetQ[0].x <= datasetP[indexesOfSplitsP[i]].x & datasetQ[0].x >= datasetP[indexesOfSplitsP[i-1]].x ){
                            rule2 = true;
                        }
                    }
                }
                
            }

            if(rule1 == false & rule2 == false){

                if(j != 0){
                    if(datasetP[indexesOfSplitsP[i]].x <= datasetQ[indexesOfSplitsQ[j]].x & datasetP[indexesOfSplitsP[i]].x >= datasetQ[indexesOfSplitsQ[j-1]].x){
                        rule4 = true;
                    }
                }
                else{
                    if(datasetP[indexesOfSplitsP[i]].x <= datasetQ[indexesOfSplitsQ[j]].x & datasetP[indexesOfSplitsP[i]].x >= datasetQ[0].x){
                        rule4 = true;
                    }
                }
                
            }

            if(rule1 == false & rule2 == false & rule4 == false){

                if(j != 0){
                    if(datasetP[indexesOfSplitsP[i]].x <= datasetQ[indexesOfSplitsQ[j-1]].x){
                        dist = abs(datasetP[indexesOfSplitsP[i]].x - datasetQ[indexesOfSplitsQ[j-1]].x);
                    
                        if( dist <= bound ){
                            rule3 = true; 
                        }
                    }
                    else{
                        if(i != 0){
                            dist = abs(datasetP[indexesOfSplitsP[i-1]].x - datasetQ[indexesOfSplitsQ[j]].x);
                    
                            if( dist <= bound ){
                                rule3 = true;  
                            }    
                        }
                        else{
                            dist = abs(datasetP[0].x - datasetQ[indexesOfSplitsQ[j]].x);
                    
                            if( dist <= bound ){
                                rule3 = true;
                            }
                        }
                    }
                    
                }
                else{
                    
                    if(datasetP[indexesOfSplitsP[i]].x <= datasetQ[0].x){
                        dist = abs(datasetP[indexesOfSplitsP[i]].x - datasetQ[0].x);
                    
                        if( dist <= bound ){
                            rule3 = true;  
                        }
                    }
                    else{
                        if(i != 0){
                            dist = abs(datasetP[indexesOfSplitsP[i-1]].x - datasetQ[indexesOfSplitsQ[j]].x);
                    
                            if( dist <= bound ){
                                rule3 = true;
                            }
    
                        }
                        else{
                            dist = abs(datasetP[0].x - datasetQ[indexesOfSplitsQ[j]].x);
                    
                            if( dist <= bound ){
                                rule3 = true; 
                            }
                            
                        }
                    }
                    
                }
                
            }

            if(readyTobreak == true){
                if(rule1 == false & rule2 == false & rule3 == false & rule4 == false){
                    readyTobreak = false;
                    break;
                }
            }
            

            if(rule1 == true || rule2 == true || rule3 == true || rule4 == true){

                neccArray.resize(count+1);
                neccArray[count].splitP = indexesOfSplitsP[i];
                neccArray[count].splitQ = indexesOfSplitsQ[j];
                
                if(i == 0){
                    neccArray[count].startP = 0;
                }
                if(j == 0){
                    neccArray[count].startQ = 0;
                }
                if(i != 0){
                    neccArray[count].startP = indexesOfSplitsP[i-1];
                }
                if(j != 0){
                    neccArray[count].startQ = indexesOfSplitsQ[j-1];
                }
                count++;
                //once found an overlap partition after that we can avoid more calculations
                readyTobreak = true;
                rule1 = false;
                rule2 = false;
                rule3 = false;
                rule4 = false;
            }
            

        }
        readyTobreak = false;
        if(i > 0){

            for(j=i-1; j >= 0; j--){
                
                //the two partitions overlap
                if(datasetQ[indexesOfSplitsQ[j]].x <= datasetP[indexesOfSplitsP[i]].x & (datasetQ[indexesOfSplitsQ[j]].x >= datasetP[indexesOfSplitsP[i-1]].x )){
                    rule1 = true;
                }
                
                if(rule1 == false){
                    if(j != 0 ){
                        //the two partitions overlap
                        if(datasetQ[indexesOfSplitsQ[j-1]].x <= datasetP[indexesOfSplitsP[i]].x & (datasetQ[indexesOfSplitsQ[j-1]].x >= datasetP[indexesOfSplitsP[i-1]].x )){
                            rule2 = true;
                        }
                    }
                    else{
                        //the two partitions overlap
                        if(datasetQ[0].x <= datasetP[indexesOfSplitsP[i]].x & (datasetQ[0].x >= datasetP[indexesOfSplitsP[i-1]].x )){
                            rule2 = true;
                        }
                    }

                }

                if(rule1 == false & rule2 == false){

                    if(j != 0){
                        if(datasetP[indexesOfSplitsP[i]].x <= datasetQ[indexesOfSplitsQ[j]].x & datasetP[indexesOfSplitsP[i]].x >= datasetQ[indexesOfSplitsQ[j-1]].x){
                            rule4 = true;
                        }
                    }
                    else{
                        if(datasetP[indexesOfSplitsP[i]].x <= datasetQ[indexesOfSplitsQ[j]].x & datasetP[indexesOfSplitsP[i]].x >= datasetQ[0].x){
                            rule4 = true;
                        }
                    }
                    
                }
                
                if(rule1 == false & rule2 == false & rule4 == false){
                    
                    if(datasetP[indexesOfSplitsP[i-1]].x >= datasetQ[indexesOfSplitsQ[j]].x){
                        dist = abs(datasetP[indexesOfSplitsP[i-1]].x - datasetQ[indexesOfSplitsQ[j]].x);
                    
                        if( dist <= bound ){
                            rule3 = true;
                        }
                    }
                    else{

                        if(j != 0){

                            dist = abs(datasetP[indexesOfSplitsP[i]].x - datasetQ[indexesOfSplitsQ[j-1]].x);
                    
                            if( dist <= bound ){
                                rule3 = true;
                            }
                        }
                        else{

                            dist = abs(datasetP[indexesOfSplitsP[i]].x - datasetQ[0].x);
                    
                            if( dist <= bound ){
                                rule3 = true;
                            }
                            
                        }
                        
                    }
                    
                }
                
                if(readyTobreak == true){
                    if(rule1 == false & rule3 == false & rule2 == false & rule4 == false){
                        readyTobreak = false;
                        break;
                    }
                }
                
                if(rule1 == true || rule3 == true || rule2 == true || rule4 == true){
                    neccArray.resize(count+1);
                    neccArray[count].splitP = indexesOfSplitsP[i];
                    neccArray[count].splitQ = indexesOfSplitsQ[j];
                    neccArray[count].startP = indexesOfSplitsP[i-1];

                    if(j == 0){
                        neccArray[count].startQ = 0;
                    }
                    else{
                        neccArray[count].startQ = indexesOfSplitsQ[j-1];
                    }
                    count++;
                    readyTobreak = true;
                    rule1 = false;
                    rule2 = false;
                    rule3 = false;
                    rule4 = false;
                }
                
            }
        }
        readyTobreak = false;      
    }   
    return neccArray;
}

void findOverlapPartitions(point_vector_t datasetP, point_vector_t datasetQ, int* indexesOfSplitsP, int* indexesOfSplitsQ, ovpairs* neccPairs, int NofSamplePartitions, int NofPartitions){

    int i;
    int j;
    int index;
    int min;
    int priority = 0; 
    bool rule1 = false;
    bool rule2 = false;
    bool rule4 = false;
    

    for(i=0; i < NofSamplePartitions; i++){

        for(j=i; j < NofSamplePartitions; j++){

            if(i == 0 ){
                
                //the two partitions overlap
                if(datasetQ[indexesOfSplitsQ[j]].x <= datasetP[indexesOfSplitsP[i]].x & datasetQ[indexesOfSplitsQ[j]].x >= datasetP[0].x){
                    rule1 = true;
                }
              
            }
            else{
                //the two partitions overlap
                if(datasetQ[indexesOfSplitsQ[j]].x <= datasetP[indexesOfSplitsP[i]].x & datasetQ[indexesOfSplitsQ[j]].x >= datasetP[indexesOfSplitsP[i-1]].x ){
                    rule1 = true;
                }

                
            }

            //the partitions overlap
            if(j != 0){
                if(i != 0){
                    if(datasetQ[indexesOfSplitsQ[j-1]].x <= datasetP[indexesOfSplitsP[i]].x & datasetQ[indexesOfSplitsQ[j-1]].x >= datasetP[indexesOfSplitsP[i-1]].x ){
                        rule2 = true;
                    }
                }
                else{
                    if(datasetQ[indexesOfSplitsQ[j-1]].x <= datasetP[indexesOfSplitsP[i]].x & datasetQ[indexesOfSplitsQ[j-1]].x >= datasetP[0].x ){
                        rule2 = true;
                    }
                }
                
            }
            else{
                if(i == 0){
                    
                    if(datasetQ[0].x <= datasetP[indexesOfSplitsP[i]].x & datasetQ[0].x >= datasetP[0].x ){
                        rule2 = true;
                    }
                }
                else{
                    
                    if(datasetQ[0].x <= datasetP[indexesOfSplitsP[i]].x & datasetQ[0].x >= datasetP[indexesOfSplitsP[i-1]].x ){
                        rule2 = true;
                    }
                }
            }
            
            if(j != 0){
                if(datasetP[indexesOfSplitsP[i]].x <= datasetQ[indexesOfSplitsQ[j]].x & datasetP[indexesOfSplitsP[i-1]].x >= datasetQ[indexesOfSplitsQ[j-1]].x){
                    rule4 = true;
                }
            }
            else{
                if(datasetP[indexesOfSplitsP[i]].x <= datasetQ[indexesOfSplitsQ[j]].x & datasetP[indexesOfSplitsP[i-1]].x >= datasetQ[0].x){
                    rule4 = true;
                }
            }
            
            if(rule1 == true){
                priority++;
            }

            if(rule2 == true){
                priority++;
            }

            if(rule4 == true){
                priority = priority + 2;
            }

            if(rule1 == false & rule2 == false & rule4 == false){
                priority = 0;
            }

            min = neccPairs[0].priority;
            index = 0;
            for(int k=1; k < NofSamplePartitions; k++){
                if(neccPairs[k].priority < min){
                    min = neccPairs[k].priority;
                    index = k;
                }
            }

            if(min == 2){
                rule1 = false;
                rule2 = false;
                rule4 = false;
                priority = 0;
                break;
            }

            if(min < priority){
                neccPairs[index].priority = priority;
                neccPairs[index].splitP = indexesOfSplitsP[i];
                neccPairs[index].splitQ = indexesOfSplitsQ[j];
                if(i == 0){
                    neccPairs[index].startP = 0;
                }
                if(j == 0){
                    neccPairs[index].startQ = 0;
                }
                if(i != 0){
                    neccPairs[index].startP = indexesOfSplitsP[i-1];
                }
                if(j != 0){
                    neccPairs[index].startQ = indexesOfSplitsQ[j-1];
                }
            }
            
            rule1 = false;
            rule2 = false;
            rule4 = false;
            priority = 0;

        }
        
        if(i > 0){

            for(j=i-1; j >= 0; j--){
                
                //the two partitions overlap
                if(datasetQ[indexesOfSplitsQ[j]].x <= datasetP[indexesOfSplitsP[i]].x & (datasetQ[indexesOfSplitsQ[j]].x >= datasetP[indexesOfSplitsP[i-1]].x )){
                    rule1 = true;
                }
                
                
                if(j != 0 ){
                    //the two partitions overlap
                    if(datasetQ[indexesOfSplitsQ[j-1]].x <= datasetP[indexesOfSplitsP[i]].x & (datasetQ[indexesOfSplitsQ[j-1]].x >= datasetP[indexesOfSplitsP[i-1]].x )){
                        rule2 = true;
                    }
                }
                else{
                    //the two partitions overlap
                    if(datasetQ[0].x <= datasetP[indexesOfSplitsP[i]].x & (datasetQ[0].x >= datasetP[indexesOfSplitsP[i-1]].x )){
                        rule2 = true;
                    }
                }

            
                if(j != 0){
                    if(datasetP[indexesOfSplitsP[i]].x <= datasetQ[indexesOfSplitsQ[j]].x & datasetP[indexesOfSplitsP[i-1]].x >= datasetQ[indexesOfSplitsQ[j-1]].x){
                        rule4 = true;
                    }
                }
                else{
                    if(datasetP[indexesOfSplitsP[i]].x <= datasetQ[indexesOfSplitsQ[j]].x & datasetP[indexesOfSplitsP[i-1]].x >= datasetQ[0].x){
                        rule4 = true;
                    }
                }
                    
                
                if(rule1 == true){
                    priority++;
                }
    
                if(rule2 == true){
                    priority++;
                }
    
                if(rule4 == true){
                    priority = priority + 2;
                }
    
                if(rule1 == false & rule2 == false & rule4 == false){
                    priority = 0;
                }
    
                min = neccPairs[0].priority;
                index = 0;
                for(int k=1; k < NofSamplePartitions; k++){
                    if(neccPairs[k].priority < min){
                        min = neccPairs[k].priority;
                        index = k;
                    }
                }

                if(min == 2){
                    rule1 = false;
                    rule2 = false;
                    rule4 = false;
                    priority = 0;
                    break;
                }

                if(min < priority){
                    neccPairs[index].priority = priority;
                    neccPairs[index].splitP = indexesOfSplitsP[i];
                    neccPairs[index].splitQ = indexesOfSplitsQ[j];
                    neccPairs[index].startP = indexesOfSplitsP[i-1];

                    if(j == 0){
                        neccPairs[index].startQ = 0;
                    }
                    else{
                        neccPairs[index].startQ = indexesOfSplitsQ[j-1];
                    }
                    
                }
                
                rule1 = false;
                rule2 = false;
                rule4 = false;
                priority = 0;
                
            }
        }
        if(min == 2){
            break;
        }      
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

