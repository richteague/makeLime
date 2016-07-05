// findvalue    - wrapper for the findcell and interpolation functions.
// findcell     - finds the bounding indices.
// linterpolate - linear interpolation.


void findcell(double c1, double c2, double c3, int *aidx, int *bidx, int *cidx,
              int *didx, int *eidx, int *fidx, int *gidx, int *hidx) {

    // Finds the appropriate cell boundaries for the given point.
    // Assumes the three coords are r = c1, z = c2 and  phi = c3.
    // The arrays are named c1arr, c2arr and c3arr respectively.
    // Returning aidx = -1 says that the point cannot be interpolated.

    double c1lower, c1upper, c3lower, c3upper;
    int i = 0;

    // Find the bounds of the phi points.

    while (c3arr[i] < c3) {
        if (i == NCELLS-1) {
            break;
        }
    }

    c3lower = c3arr[i];
    c3higher = c3arr[0];

    // For each of the phi slices, do as in the 2D case.
    // Find the bounds of the radial points, lower phi position.

    i = 0;
    for (i=0; i<(NCELLS-1); i=i+1) {
        if (c3arr[i] == c3lower && c3arr[i-1] == c3lower) {
            if ((c1arr[i]-c1) * (c1arr[i-1]-c1) < 0.) {
                break;
            }
        }
    }

    c1upper = c1arr[i];
    c1lower = c1arr[i-1];

   // Find the bounds of the vertical points at the lower radial position, lower phi position.

    i = 0;
    for(i=0; i<(NCELLS-1); i=i+1) {
        if (c3arr[i] == c3lower && c1arr[i] == c1lower && c1arr[i-1] == c1lower) {
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

    // Find the bounds of the vertical points at the upper radial position, lower phi position.

    i = 0;
    for(i=0; i<(NCELLS-1); i=i+1) {
        if (c3arr[i] == c3lower && c1arr[i] == c1upper && c1arr[i-1] == c1upper) {
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
    *cidx = i-1;

    // Find the bounds of the radial points, upper phi position.

    i = 0
    for (i=0; i<(NCELLS-1); i=i+1) {
        if (c3arr[i] == c3upper && c3arr[i-1] == c3upper) {
            if ((c1arr[i]-c1) * (c1arr[i-1]-c1) < 0.) {
                break;
            }
        }
    }

    c1upper = c1arr[i];
    c1lower = c1arr[i-1];

    // Find the bounds of the vertical points at the lower radial position, lower phi position.

    i = 0;
    for(i=0; i<(NCELLS-1); i=i+1) {
        if (c3arr[i] == c3upper && c1arr[i] == c1lower && c1arr[i-1] == c1lower) {
            if ((c2arr[i]-c2)*(c2arr[i-1]-c2) < 0.) {
                break;
            }
        }
    }

    if (c1arr[i] != c1lower) {
        *aidx = -1;
        return;
    }

    *fidx = i;
    *eidx = i-1;

    // Find the bounds of the vertical points at the upper radial position, lower phi position.

    i = 0;
    for(i=0; i<(NCELLS-1); i=i+1) {
        if (c3arr[i] == c3upper && c1arr[i] == c1upper && c1arr[i-1] == c1upper) {
            if ((c2arr[i]-c2)*(c2arr[i-1]-c2) < 0.) {
                break;
            }
        }
    }

    if (c1arr[i] != c1upper) {
        *aidx = -1;
        return;
    }

    *hidx = i;
    *gidx = i-1;

}


double linterpolate(double c1, double c2, double c3, int aidx, int bidx, int cidx, int didx, 
                    int eidx, int fidx, int gidx, int hidx, const double arr[NCELLS]) {

    // Three dimensional, linear interpolation.

    double A, B, C, D, f;

    // Lower phi component.

    f = (c2 - c2arr[aidx]) / (c2arr[bidx] - c2arr[aidx]);
    A = arr[aidx] * (1. - f) + arr[bidx] * f;

    f = (c2 - c2arr[cidx]) / (c2arr[didx] - c2arr[cidx]);
    B = arr[cidx] * (1. - f) + arr[didx] * f;

    f = (c1 - c1arr[aidx]) / (c1arr[cidx] - c1arr[aidx]);
    C = A * (1. - f) + B * f;

    // Upper phi component.

    f = (c2 - c2arr[eidx]) / (c2arr[fidx] - c2arr[eidx]);
    A = arr[eidx] * (1. - f) + arr[fidx] * f;

    f = (c2 - c2arr[gidx]) / (c2arr[hidx] - c2arr[gidx]);
    B = arr[gidx] * (1. - f) + arr[hidx] * f;

    f = (c1 - c1arr[eidx]) / (c1arr[gidx] - c1arr[eidx]);
    D = A * (1. - f) + B * f;

    // Interpolate phi components.

    f = (c3 - c3arr[aidx]) / (c3arr[eidx] - c3arr[aidx]);
    return C * (1. - f) + D * f

}


double findvalue(double c1, double c2, double c3, const double arr[NCELLS]){

    // Finds the bounding cells and linerally interpolates their value.

    int aidx, bidx, cidx, didx, eidx, fidx, gidx, hidx;

    findcell(c1, c2, c3, &aidx, &bidx, &cidx, &didx, &eidx, &fidx, &gidx, &hidx);

    if (aidx >= 0) {
        return linterpolate(c1, c2, c3, aidx, bidx, cidx, didx,
                            eidx, fidx, gidx, hidx, arr);
    } else {
        return -1.;
    }
}
