// findvalue    - wrapper for the findcell and interpolation functions.
// findcell     - finds the bounding indices.
// linterpolate - linear interpolation.


void findcell(double c1, double c2, int *aidx, int *bidx, int *cidx, int *didx){

    // In this function, returning *aidx < 0 means that no cell has been found.

    double c1lower, c1upper;
    int i = 0;

    // Find the bounding radial points such that c1lower < r <= c1upper.

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

    /*
        Find the bounds of the vertical points at the lower radial position.
        This is found by (z1 - z) / (z0 - z) < 0. At the end of the loop, if
        c1arr[i] != r then we have overshot. In case the point is successfully
        found in the upper radial position then we continue and leave
        aidx = bidx.
    */

    i = 0;
    for (i=0; i<(NCELLS-1); i=i+1) {
        if (c1arr[i] == c1lower && c1arr[i-1] == c1lower) {
            if ((c2arr[i]-c2)*(c2arr[i-1]-c2) < 0.) {
                *bidx = i;
                *aidx = i-1;
                break;
            }
        }
        if (c1arr[i] > c1lower) {
            *bidx = i-1;
            *aidx = i-1;
            break;
        }
    }

    /*
        At this point we do not know if z is within the bounds of the first
        vertical column (aidx != bidx) or if it is above (aidx == bidx).
    */

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
    *cidx = i-1;

    /*
        At this point we know that z is in the bounds of the upper vertical
        column, but we do not know if it is in the bounds of the lower vertical
        column (aidx != bidx if so). If it is also in, then we can return the
        four indices as usual. If not, we calculate to make sure the projected
        value is
    */

    if (*aidx == *bidx) {
        if (c2 > c1arr[*didx] * c2arr[*aidx] / c1arr[*aidx]){
            *aidx = -1;
            return;
        }
    }

    /*
        If we get to here, then interpolation must be possible.
    */

    return;

}


double linterpolate(double c1, double c2, int aidx, int bidx, int cidx,
                    int didx, const double arr[NCELLS]) {

    /*
        Two dimensional, linear interpolation. Consider two cases: bi-linear
        if aidx != bidx, or radial otherwise.
    */

    double A, B, f;

    // Standard bi-linear interpolation.

    if (aidx != bidx) {
        f = (c2 - c2arr[aidx]) / (c2arr[bidx] - c2arr[aidx]);
        A = arr[aidx] * (1. - f) + arr[bidx] * f;
        f = (c2 - c2arr[cidx]) / (c2arr[didx] - c2arr[cidx]);
        B = arr[cidx] * (1. - f) + arr[didx] * f;
        f = (c1 - c1arr[aidx]) / (c1arr[cidx] - c1arr[aidx]);
        return A * (1. - f) + B * f;
    }

    // Overshot lower radial position.

    double c3;

    c3 = c2arr[aidx] * c1arr[didx] / c1arr[aidx];
    f = (c3 - c2arr[cidx]) / (c2arr[didx] - c2arr[cidx]);
    B = arr[cidx] * (1. - f) + arr[didx] * f;
    f = (c1 - c1arr[aidx]) / (c1arr[cidx] - c1arr[aidx]);
    return arr[aidx] * (1. - f) + B * f;

}


double findvalue(double c1, double c2, double c3, const double arr[NCELLS]){

    /*
        Finds the bounding cells and linearlly interpolates their value.
        If aidx < 0 then we cannot interpolate this point in the provided model
        grid.
    */

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
