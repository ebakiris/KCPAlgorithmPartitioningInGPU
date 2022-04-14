void BruteForceAlgorithm(point* datasetP, point* datasetQ, int NofCP, int nofPointsP, int nofPointsQ){


    bool flag;
    int i;
    int j;
    int index;
    double dist;
    closestpairs* localMaxKHeap;

    localMaxKHeap = (closestpairs*)malloc(NofCP*sizeof(closestpairs));

    //INITIALIZATION
    for(i = 0; i < NofCP; i++){
        localMaxKHeap[i].dist = DBL_MIN;
    }

    for(i = 0; i < nofPointsP; i++){

        for(j = 0; j < nofPointsQ; j++){
            
            dist = EuclideanHost(datasetP[i].x, datasetP[i].y, datasetQ[j].x, datasetQ[j].y);

            flag = isFullHost(localMaxKHeap, NofCP);
            index = insertDistanceHost(dist, localMaxKHeap, flag, NofCP);
            
            if(index != -1){
                localMaxKHeap[index].p = datasetP[i];
                localMaxKHeap[index].q = datasetQ[j];
            }
        }
    }

    for(i = 0; i < NofCP; i++){

        printf("The pair is: p (%f, %f) and q (%f, %f) with Distance: %.10lf\n", localMaxKHeap[i].p.x, localMaxKHeap[i].p.y, localMaxKHeap[i].q.x, localMaxKHeap[i].q.y, localMaxKHeap[i].dist);
        
    }

}