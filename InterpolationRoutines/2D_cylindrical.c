// findvalue    - wrapper for the findcell and interpolation functions.
// findcell     - finds the bounding indices.
// linterpolate - linear interpolation.


void findcell(double c1, double c2, int *aidx, int *bidx, int *cidx, int *didx){

    double c1lower, c1upper;
    int i = 0;

    // Find the bounding radial points.

    i = 0;
    while (c1arr[i] < c1) {
        if (i == NCELLS - 1) {
            *aidx = -1;
            return;
        }
        i++;
    }    

    c1lower = c1arr[i-1];
    c1upper = c1arr[i];

    // Find the bounds of the vertical points at the lower radial position.

    i = 0;
    for (i=0; i<(NCELLS-1); i=i+1) {
        if (c1arr[i] == c1lower && c1arr[i-1] == c1lower) {
            if ((c2arr[i]-c2)*(c2arr[i-1]-c2) < 0.) {
                break;
            }
        }
    }

    if (c1arr[i] != c1lower) {
        *aidx = -1;
        return;
    }

    *bidx = i;
    *aidx = i-1;

    // Find the bounds of the vertical points at the upper radial position.

    i = 0;
    for (i=0; i<(NCELLS-1); i=i+1) {
        if (c1arr[i] == c1upper && c1arr[i-1] == c1upper) {
            if ((c2arr[i]-c2)*(c2arr[i-1]-c2) < 0.) {
                break;
            }
        }
    }

    if (c1arr[i] != c1upper) {
        *aidx = -1;
        return;
    }

    *didx = i;
    *cidx = i -1;

}


double linterpolate(double c1, double c2, int aidx, int bidx, int cidx, 
                    int didx, const double arr[NCELLS]) {

    // Two dimensional, linear interpolation.

    double A, B, f;

    f = (c2 - c2arr[aidx]) / (c2arr[bidx] - c2arr[aidx]);
    A = arr[aidx] * (1. - f) + arr[bidx] * f;

    f = (c2 - c2arr[cidx]) / (c2arr[didx] - c2arr[cidx]);
    B = arr[cidx] * (1. - f) + arr[didx] * f;

    f = (c1 - c1arr[aidx]) / (c1arr[cidx] - c1arr[aidx]);
    return A * (1. - f) + B * f;

}


double findvalue(double c1, double c2, double c3, const double arr[NCELLS]){

    // Finds the bounding cells and linearlly interpolates their value.

    double value;
    int aidx, bidx, cidx, didx;

    findcell(c1, c2, &aidx, &bidx, &cidx, &didx);

    if (aidx >= 0) {
        value = linterpolate(c1, c2, aidx, bidx, cidx, didx, arr);
        return value;
    } else {
        return -1.;
    }

}
