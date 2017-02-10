// findvalue    - wrapper for the findcell and interpolation functions.
// findcell     - finds the bounding indices.
// linterpolate - linear interpolation.

double projected_value(double c1, double c2, int aidx, int cidx){

    /*
        Calculate the projection of the z value on the outer column.
    */

    if (c1 == c1arr[aidx]) {
        return -1.0;
    }

    double z_proj;
    z_proj = c1arr[cidx] - c1arr[aidx];
    z_proj *= c2 - c2arr[aidx];
    z_proj /= c1 - c1arr[aidx];
    return z_proj + c2arr[aidx];
}

void verticalbounds(double c1val, double c2, int *lidx, int *uidx){
    
    /*
        Find the bounding cells, lidx and uidx, of 'c2' in a vertical column
        at radius c1val. Returns lidx == uidx if not found in column. 
    */

    int i = 0;
    for (i=0; i<(NCELLS-1); i=i+1) {
        if (c1arr[i] == c1val && c1arr[i-1] == c1val) {
            if ((c2arr[i]-c2)*(c2arr[i-1]-c2) < 0.) {
                *uidx = i;
                *lidx = i-1;
                return;
            }
        } else if (c1arr[i] > c1val) {
            *uidx = i-1;
            *lidx = i-1;
            break;
        }
    }
    return;
}

void findcell(double c1, double c2, int *aidx, int *bidx, int *cidx, int *didx){

    // In this function, returning *aidx < 0 means that no cell has been found.

    double c1lower, c1upper;
    double ylim, m, c;
    int atemp, btemp, ctemp, dtemp, zmax;
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
        c1arr[i] != r then we have overshot. In the case the vertical point
        is not found, then aidx = bidx or cidx = didx.
    */

    atemp = -1;
    verticalbounds(c1lower, c2, &atemp, &btemp);
    verticalbounds(c1upper, c2, &ctemp, &dtemp);

    *aidx = atemp;
    *bidx = btemp;
    *cidx = ctemp;
    *didx = dtemp;

    if (*cidx == *didx) {

        // No interpolation possible.
        *aidx = -1;
        return;

    } else if (*aidx != *bidx) {

        // Bilinear interpolation possible.
        return;

    } else {

        // Check if 'tripolation' is possible.
        for (i=0; i<NCELLS-1; i=i+1){
            if (c1arr[i] > c2){
                break;
            }
        }

        if (c1arr[i] == c2) {
            zmax = i;
        } else if (c1arr[i] > c2) {
            zmax = i - 1;
        } else {
            *aidx = -1;
            return;
        }

        m = c2arr[zmax] - c2arr[*aidx];
        m /= c1arr[zmax] - c1arr[*aidx];
        c = c2arr[zmax] - m * c1arr[zmax];

        if (c2 > (m * c1 + c)) {
            *aidx = -1;
            return;
        }
    
    }
    
    // Find the projected z value on the outer column
    // and return those bounding indicies.
    
    // double zproj;
    //zproj = projected_value(c1, c2, atemp, ctemp);
    //    verticalbounds(c1upper, zproj, &ctemp, &dtemp);
    
    *aidx = -1;
    return;

}


double linterpolate(double c1, double c2, int aidx, int bidx, int cidx,
                    int didx, const double arr[NCELLS]) {

    /*
        Two dimensional, bilinear interpolation.
    */

    double A, B, f;
    f = (c2 - c2arr[aidx]) / (c2arr[bidx] - c2arr[aidx]);
    A = arr[aidx] * (1. - f) + arr[bidx] * f;
    f = (c2 - c2arr[cidx]) / (c2arr[didx] - c2arr[cidx]);
    B = arr[cidx] * (1. - f) + arr[didx] * f;
    f = (c1 - c1arr[aidx]) / (c1arr[cidx] - c1arr[aidx]);
    return A * (1. - f) + B * f;

}



double tripolate(double c1, double c2, int aidx, int cidx, int didx, const double arr[NCELLS]) {

    /*
        Interpolating the triangular region above the maximum of the inner
        radial grid, but still within the bounds of the outer grid. aidx
        specifies the top-most point of the c1lower column, while 
    */

    double A, B, f;
    double z_proj;

    z_proj = projected_value(c1, c2, aidx, cidx);
    A = arr[aidx];
    f = (c2 - c2arr[cidx]) / (c2arr[didx] - c2arr[cidx]);
    B = arr[cidx] * (1. - f) + arr[didx] * f;
    f = (c1 - c1arr[aidx]) / (c1arr[cidx] - c1arr[aidx]);
    return A * (1. - f) + B * f;

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
        if (value < 0){
            return -1.;
        }
        return value;
    } else {
        return -1.;
    }

}
