void findcell(double rad, double alt, int *aidx, 
              int *bidx, int *cidx, int *didx){

    double rlower, rupper;
    int i = 0;
   
    // Find the bounds of the radial points.
    while (rvals[i] < rad) {
        if (i == NCELLS-1) {
            *aidx = -1;
            return;
        }
        i++;
    }

    rlower = rvals[i-1];
    rupper = rvals[i];
    
    // Find the bounds of the vertical points at the lower radial position.
    i = 0;
    for(i=0; i<(NCELLS-1); i=i+1) {
        if (rvals[i] == rlower && rvals[i-1] == rlower) {
            if ((zvals[i]-alt)*(zvals[i-1]-alt) < 0.) {
                break;
            }
        }
    }
    
    if (rvals[i] != rlower) {
        *aidx = -1;
        return;
    }

    *bidx = i;
    *aidx = i-1; 

    // Find the bounds of the vertical points at the upper radial position.
    i = 0;
    for(i=0; i<(NCELLS-1); i=i+1) {    
        if (rvals[i] == rupper && rvals[i-1] == rupper) {
            if ((zvals[i]-alt)*(zvals[i-1]-alt) < 0.) {
                break;
            }
        }
    }  
    
    if (rvals[i] != rupper) {
        *aidx = -1;
        return;
    }  
      
    *didx = i;
    *cidx = i-1;         
}
