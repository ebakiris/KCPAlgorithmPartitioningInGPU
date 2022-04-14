#include <cmath>

__device__ bool isFull(void* localMaxKHeap, int NofCP, bool phase2){

    int i;
    bool value = true;

    if(phase2 == false){
        double* local = (double*) localMaxKHeap;
        for(i = 0; i < NofCP; i++){    
            
            if(local[i] <= DBL_MIN){
                value = false;
                break;
            }
        }  
    }
    else{
        closestpairs* local = (closestpairs*) localMaxKHeap;
        for(i = 0; i < NofCP; i++){    
            if(local[i].dist <= DBL_MIN){
                value = false;
                break;
            }
        }      
    }
    
    return value;

}

__device__ double Euclidean(double x1, double y1, double x2, double y2)
{
	double x = x1 - x2; //calculating number to square in next step
	double y = y1 - y2;
	double dist;

	dist = pow(x, 2) + pow(y, 2);       //calculating Euclidean distance
	dist = sqrt(dist);                  

	return dist;
}


__device__ double maxDelta(void* localMaxKHeap, int NofCP, bool phase2){

    if(!phase2){
        double* local = (double*) localMaxKHeap;
        double max = local[0];

        for(int i=1; i < NofCP; i++){

            if(local[i] > max){
                max = local[i];
            }
        }
        return max;
    }
    else{
        closestpairs* local = (closestpairs*) localMaxKHeap;
        double max = local[0].dist;

        for(int i=1; i < NofCP; i++){

            if(local[i].dist > max){
                max = local[i].dist;
            }
        }
        return max;
    }
}

__device__ int insertDistance(double dist, void* localMaxKHeap, bool isfull, int NofCP, bool phase2){

    double max;
    int i;
    
    
    if(isfull == false){

        if(phase2 == false){
            double* local = (double*) localMaxKHeap;
            i = 0;
            while(!(local[i] <= DBL_MIN)){
                i++;
            }
            local[i] = dist;
            return i;
        }
        else{
            closestpairs* local = (closestpairs*) localMaxKHeap;
            i = 0;
            while(!(local[i].dist <= DBL_MIN)){
                i++;
            }
            local[i].dist = dist;
            return i;
        }
    }
    else{
        
        if(phase2 == false){

            double* local = (double*) localMaxKHeap;
            
            int index = 0;
            max = local[0];

            for(i=1; i < NofCP; i++){

                if(local[i] > max){
                    max = local[i];
                    index = i;
                }
            }

            if(dist < max & dist > DBL_MIN){

                local[index] = dist;
                return index;
                
            }
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

            if(dist < max & dist > DBL_MIN){

                local[index].dist = dist;
                return index;
                
            }
        }
        

    }
    return -1;
}
