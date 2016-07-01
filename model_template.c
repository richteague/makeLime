/*
 *  model.c
 *  LIME, The versatile 3D line modeling tool 
 *
 *  Created by Christian Brinch on 11/05/07.
 *  Copyright 2006-2013, Christian Brinch, 
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "lime.h"
#include "math.h"
#include "stdio.h"
#include <stdlib.h>

// Include

// Chemical Model

// Stellar Mass

// Mach

// Grid Cells

// Ortho Ratio

void input(inputPars *par, image *img)
{
// Input Parameters

// Images
// 0

}



/*
                -- Interpolation Functions -- 
                
    findcell     - finds the bounding cell indicies assuming a 1+1 grid.
    linterpolate - linearally interpolates between the boudning cells.
    findvalue    - wrapper for the two functions.

*/

void findcell(double rad, double alt, int *aidx, 
              int *bidx, int *cidx, int *didx)
{

    double rlower, rupper;
    int i = 0;
   
    // Find the bounds of the radial points.
    while (rvals[i] < rad)
    {
        if (i == NCELLS)
        {
            *aidx = -1;
            return;
        }
        i++;
    }
    rlower = rvals[i-1];
    rupper = rvals[i];
    
    // Find the bounds of the vertical points at the lower radial position.
    i = 0;
    for(i=0; i<(NCELLS-1); i=i+1)
    {    
        if (rvals[i] == rlower && rvals[i-1] == rlower)
        {
            if ((zvals[i]-alt)*(zvals[i-1]-alt) < 0.)
            {
                break;
            }
        }
    }
    
    if (rvals[i] != rlower)
    {
        *aidx = -1;
        return;
    }
    *bidx = i;
    *aidx = i-1; 

    // Find the bounds of the vertical points at the upper radial position.
    i = 0;
    for(i=0; i<(NCELLS-1); i=i+1)
    {    
        if (rvals[i] == rupper && rvals[i-1] == rupper)
        {
            if ((zvals[i]-alt)*(zvals[i-1]-alt) < 0.)
            {
                break;
            }
        }
    }  
    
    if (rvals[i] != rupper)
    {
        *aidx = -1;
        return;
    }  
      
    *didx = i;
    *cidx = i-1;         
}

double linterpolate(double rad, double alt, int aidx, int bidx,
                    int cidx, int didx, const double arr[NCELLS])
{	

    double A, B, f;
    f = (alt - zvals[aidx]) / (zvals[bidx] - zvals[aidx]);
    A = arr[aidx] * (1. - f) + arr[bidx] * f;
    f = (alt - zvals[cidx]) / (zvals[didx] - zvals[cidx]);
    B = arr[cidx] * (1. - f) + arr[didx] * f;
    f = (rad - rvals[aidx]) / (rvals[cidx] - rvals[aidx]);
    return A * (1. - f) + B * f;

}

double findvalue(double rad, double alt, const double arr[NCELLS])
{

    double value;
    int aidx, bidx, cidx, didx;
    
    findcell(rad, alt, &aidx, &bidx, &cidx, &didx);
    if (aidx >= 0)
    {
        value =  linterpolate(rad, alt, aidx, bidx, cidx, didx, arr);
    } else 
    {
        value = -1.;
    }
    
    return value;
}



/*
               -- Model Functions --
               
    density     - returns the density of colliders (H2).
    temperature - returns the gas and dust temperatures.
    abundance   - returns the molecular abundance relative to collider density.
    doppler     - returns the Doppler broadening parameter.
    gassIIdust  - returns the gas to dust ratio.
    velocity    - returns the three cartestian velocity components.
    
*/

void density(double x, double y, double z, double *density){	
    // Density
}

void temperature(double x, double y, double z, double *temperature){
    // Temperature
}

void abundance(double x, double y, double z, double *abundance)
{
    // Abundance
}

void doppler(double x, double y, double z, double *doppler)
{
    // Turbulence
}

void gasIIdust(double x, double y, double z, double *gtd)
{
    // Gas To Dust
}

void velocity(double x, double y, double z, double *velocity)
{
    // Velocity
}
