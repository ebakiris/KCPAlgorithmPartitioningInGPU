//N is the size of the partition P, and M the size of the partition Q
//beginP is the index of the beginning point of partition P
//beginQ is the index of the beginning point of partition Q


void ReverseRunPlaneSweepSerial(point* datasetP, point* datasetQ, int nofPointsP, int nofPointsQ, int NofCP, double* globalMaxKHeap)
{
    
    int N ;
    int M ;
    int beginP ;
    int beginQ ;
    int i ;
    int j ;
    bool cont ;
    int k;
    double dx;
    double d;
    double delta;
    double* localMaxKHeap;
    bool flag ;

    localMaxKHeap = (double*)malloc(NofCP*sizeof(double));

    for(int l=0; l < NofCP; l++){
        localMaxKHeap[l] = DBL_MIN;
    }

    
    N = nofPointsP;
    M = nofPointsQ;
    beginP = 0;
    beginQ = 0;
    

    i = beginP;
    j = beginQ;
    cont = true;
    
    //sentinel points for simpler stoping conditions
    double sentinelP = DBL_MAX;
    double sentinelQ = DBL_MAX;
    bool flag2 = true;
    bool flag3 = true;

    int index;
    
    delta = 0.05;
    
    if(datasetP[N-1].x <= datasetQ[beginQ].x ){
        i = N;  //the sets do not overlap
    }
  
    if(datasetQ[M-1].x <= datasetP[beginP].x ){
        j = M;  //the sets do not overlap
    }

    int leftp = -1;
    int leftq = -1; //comparisons start at P[leftp+1], Q[leftq+1]

    //Main algorithm P[i](Q[j]):Start of next P-run(Q-run)
    while(cont){

        if(i == N & j == M){
                    
            flag2 = false;
            
        }
        if(i != N & j == M){
          
            flag2 = datasetP[i].x < sentinelQ;
            
        }
        if(i == N & j != M){
            
            flag2 = sentinelP < datasetQ[j].x ;
       
        }
        if(i != N & j != M){

            flag2 = datasetP[i].x < datasetQ[j].x;
            
        }
        
        if(flag2) {  //the active run is from the P set
            while( flag2 ) {   //while active run unfinished. P[i]: ref_point
                if(j - 1 == leftq) { //Q[j-1]: last cur_point - rule 3
                    i++;
                    break;
                }
                for(k = j-1; k >= leftq + 1;k--){    //Q[k]: cur_point
                    
                   
                    flag = isFullHost(localMaxKHeap, NofCP);
                        
                    if ( flag == false ){
                        //calculate distance between P[i] and Q[k]
                        d = EuclideanHost(datasetP[i].x, datasetP[i].y, datasetQ[k].x, datasetQ[k].y);
                        //insert distance into localMaxKHeap
                        insertDistanceHost(d, &localMaxKHeap[0], false, NofCP);
                    }
                    
                    if(flag == true){
                        //calculate x-distance between P[i] and Q[k]
                        dx = abs(datasetQ[k].x - datasetP[i].x) ;
                        
                        if(dx >= delta ) {     //delta is a value depending the version of rrps used(sliding window,semi-circle,etc)
                            leftq = k;
                            break;
                        }
                        //calculate distance d between ref-point(P[i]) and cur_point(Q[k])
                        d = EuclideanHost(datasetP[i].x, datasetP[i].y, datasetQ[k].x, datasetQ[k].y);
                        
                        if(d < delta){
                            index = insertDistanceHost(d, &localMaxKHeap[0], true, NofCP);
                            if(index != -1){
                                delta = maxDelta(localMaxKHeap, NofCP);
                            }
                        }
                    }
                }
                i++;
                if(i == N & j == M){
                    
                    flag2 = false;
                    
                }
                if(i != N & j == M){
                  
                    flag2 = datasetP[i].x < sentinelQ;
                    
                }
                if(i == N & j != M){
                    
                    flag2 = sentinelP < datasetQ[j].x ;
               
                }
                if(i != N & j != M){
        
                    flag2 = datasetP[i].x < datasetQ[j].x;
                    
                }
            }
        }     
        else if(j < M){     //the active run is from the (unfinished) Q set 
            
            sentinelP = datasetQ[M-1].x + 1;        //P[N] should be < Q[M], since...
            if(i == N){
                
                flag3 = datasetQ[j].x <= sentinelP ;
                
            }
            else{
                
                flag3 = datasetQ[j].x <= datasetP[i].x;
                
            }
            
            while(flag3){  //while active run unfinished. Q[j]: ref_point
                if(i-1 == leftp){                    //P[i-1]: last cur_point 
                    j++;
                    break;
                }
                for(k=i-1; k >= leftp+1; k--){     //P[k]: cur_point
                    
                    flag = isFullHost(localMaxKHeap, NofCP);
                    
                    if ( flag == false ){

                        //calculate distance d between ref_point (Q[j]) and cur_point(P[k])
                        
                        d = EuclideanHost(datasetQ[j].x, datasetQ[j].y, datasetP[k].x, datasetP[k].y);
                        //insert distance into localMaxKHeap
                        insertDistanceHost(d, &localMaxKHeap[0], false, NofCP);
                        
                    }
                    if( flag == true ){
                        //calculate x-distance dx between ref_point (Q[j]) and cur_point(P[k])
                        dx = abs(datasetP[k].x - datasetQ[j].x) ;
                        
                        if(dx >= delta ){
                            leftp = k;
                            break;
                        }
                        //calculate distance d between ref_point (Q[j]) and cur_point(P[k])
                        d = EuclideanHost(datasetQ[j].x, datasetQ[j].y, datasetP[k].x, datasetP[k].y);
                        if(d < delta){
                            //insert distance into localMaxKHeap
                            index = insertDistanceHost(d, &localMaxKHeap[0], true, NofCP);
                            if(index != -1){
                                delta = maxDelta(localMaxKHeap, NofCP);
                            }
                        }
                    }
                }
                j++;    //update ref-point Q[j]
                if(j == M & i == N){
                   
                    flag3 = false;
                    
                }
                if(j == M & i != N){
                    
                    flag3 = sentinelQ <= datasetP[i].x;

                }
                if(j != M & i != N){
                    flag3 = datasetQ[j].x <= datasetP[i].x;
                    
                }
                
                if(j != M & i == N){
                    flag3 = datasetQ[j].x <= sentinelP ;
                   
                }
                
            }
            sentinelP = sentinelQ;          //revert the P sentinel at the maximum real X-value
        }
        else{
            cont = false;
        }
    } 

    
    for(int i=0; i < NofCP; i++){
        globalMaxKHeap[i] = localMaxKHeap[i];
    }
    
}