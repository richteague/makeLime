// findvalue    - wrapper for the findcell and interpolation functions.
// findcell     - finds the bounding indices.
// linterpolate - linear interpolation.

double findvalue(double c1, double c2, double c3, const double arr[NCELLS]){

    // Finds the bounding cells and linearlly interpolates them.

    double value;
    int aidx, bidx, cidx, didx;

    findcell(cone, ctwo, &aidx, &bidx, &cidx, &didx);
    if (aidx >= 0) {
        value = linterpolate(cone, ctwo, aidx, bidx, cidx, didx, arr);
    } else {
        value = -1.;
    }

    return value;

}

