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

// Chemical Model

// Stellar Mass

// Grid Cells

void
input(inputPars *par, image *img){
// Input Parameters

// Images
// 0

}

/******************************************************************************/
void
findcell(double rad, double alt, int *aidx, int *bidx, int *cidx, int *didx){
    
    // Find the indicies of the bounding cell with which to interpolate.
    // If the value is not in the grid, return aidx = -1, which tells the 
    // interpolation to return the default value.

    double rlower, rupper;
    int i = 0;
   
    // Find the bounds of the radial points.
    while (rvals[i] < rad){
        if (i == NCELLS){
            *aidx = -1;
            return;
        }
        i++;
    }
    
    rlower = rvals[i-1];
    rupper = rvals[i];
    
    // Find the bounds of the vertical points at the lower radial position.
    i = 0;
    for(i=0; i<(NCELLS-1); i=i+1){    
        if (rvals[i] == rlower && rvals[i-1] == rlower){
            if ((zvals[i]-alt)*(zvals[i-1]-alt) < 0.){
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
    for(i=0; i<(NCELLS-1); i=i+1){    
        if (rvals[i] == rupper && rvals[i-1] == rupper){
            if ((zvals[i]-alt)*(zvals[i-1]-alt) < 0.){
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

/******************************************************************************/
void
density(double x, double y, double z, double *density){	
 
    // Declare variables here.
    double rad = sqrt(x*x + y*y)/AU;
    double alt = fabs(z)/AU;
    double A, B, f;
    int aidx, bidx, cidx, didx;
    
    // Default value.
    density[0] = 1e-30;

    // This returns the indicies of the bounding cells.
    findcell(rad, alt, &aidx, &bidx, &cidx, &didx);

    // If aidx >= 0 then we can interpolate. Otherwise return the default value.
    if (aidx >= 0){
    
        // Linearly interpolate the values here.
        f = (alt - zvals[aidx]) / (zvals[bidx] - zvals[aidx]);
        A = dens[aidx] * (1. - f) + dens[bidx] * f;
        f = (alt - zvals[cidx]) / (zvals[didx] - zvals[cidx]);
        B = dens[cidx] * (1. - f) + dens[didx] * f;
        f = (rad - rvals[aidx]) / (rvals[cidx] - rvals[aidx]);
        density[0] = A * (1. - f) + B * f;
        
        // Check to make sure this isn't below the default value.
        if (density[0] < 1e-30){
	        density[0] = 1e-30;
        }
    }    
}

/******************************************************************************/

void
temperature(double x, double y, double z, double *temperature){

    // Declare variables here.
    double rad = sqrt(x*x + y*y)/AU;
    double alt = fabs(z)/AU;
    double A, B, f;
    int aidx, bidx, cidx, didx;
    
    // Default value.
    // Note: temperature[0] is gas temp and temperature[1] is dust temp.
    temperature[0] = 2.7;
    temperature[1] = 2.7;

    // This returns the indicies of the bounding cells.
    findcell(rad, alt, &aidx, &bidx, &cidx, &didx);

    // If aidx >= 0 then we can interpolate. Otherwise return the default value.
    if (aidx >= 0){
    
        // Linearly interpolate the values here.
        f = (alt - zvals[aidx]) / (zvals[bidx] - zvals[aidx]);
        A = temp[aidx] * (1. - f) + temp[bidx] * f;
        f = (alt - zvals[cidx]) / (zvals[didx] - zvals[cidx]);
        B = temp[cidx] * (1. - f) + temp[didx] * f;
        f = (rad - rvals[aidx]) / (rvals[cidx] - rvals[aidx]);
        temperature[0] = A * (1. - f) + B * f;
        
        // Check to make sure this isn't below the default value.
        if (temperature[0] < 2.7){
	        temperature[0] = 2.7;
        }
        
        // Linearly interpolate the values here.
        f = (alt - zvals[aidx]) / (zvals[bidx] - zvals[aidx]);
        A = dtemp[aidx] * (1. - f) + dtemp[bidx] * f;
        f = (alt - zvals[cidx]) / (zvals[didx] - zvals[cidx]);
        B = dtemp[cidx] * (1. - f) + dtemp[didx] * f;
        f = (rad - rvals[aidx]) / (rvals[cidx] - rvals[aidx]);
        temperature[1] = A * (1. - f) + B * f;
        
        // Check to make sure this isn't below the default value.
        if (temperature[1] < 2.7){
	        temperature[1] = 2.7;
        }

    }

}

/******************************************************************************/

void
abundance(double x, double y, double z, double *abundance){

    // Declare variables here.
    double rad = sqrt(x*x + y*y)/AU;
    double alt = fabs(z)/AU;
    double A, B, f;
    int aidx, bidx, cidx, didx;
    
    // Default value.
    // Relative abundance to the main collision partner described by density[0].
    abundance[0] = 0.0;

    // This returns the indicies of the bounding cells.
    findcell(rad, alt, &aidx, &bidx, &cidx, &didx);

    // If aidx >= 0 then we can interpolate. Otherwise return the default value.
    if (aidx >= 0){
    
        // Linearly interpolate the values here.
        f = (alt - zvals[aidx]) / (zvals[bidx] - zvals[aidx]);
        A = abund[aidx] * (1. - f) + abund[bidx] * f;
        f = (alt - zvals[cidx]) / (zvals[didx] - zvals[cidx]);
        B = abund[cidx] * (1. - f) + abund[didx] * f;
        f = (rad - rvals[aidx]) / (rvals[cidx] - rvals[aidx]);
        abundance[0] = A * (1. - f) + B * f;
        
        // Check to make sure this isn't below the default value.
        if (abundance[0] < 0.0){
	        abundance[0] = 0.0;
        }

    }
}

/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){

	// Assume a constant value.
    *doppler = 50.;

}

/******************************************************************************/

void
gasIIdust(double x, double y, double z, double *gtd){

    // Declare variables here.
    double rad = sqrt(x*x + y*y)/AU;
    double alt = fabs(z)/AU;
    double A, B, f;
    int aidx, bidx, cidx, didx;
    
    // Default value.
    // Gas to dust ratio.
    *gtd = 100.;

    // This returns the indicies of the bounding cells.
    findcell(rad, alt, &aidx, &bidx, &cidx, &didx);

    // If aidx >= 0 then we can interpolate. Otherwise return the default value.
    if (aidx >= 0){
    
        // Linearly interpolate the values here.
        f = (alt - zvals[aidx]) / (zvals[bidx] - zvals[aidx]);
        A = gastodust[aidx] * (1. - f) + gastodust[bidx] * f;
        f = (alt - zvals[cidx]) / (zvals[didx] - zvals[cidx]);
        B = gastodust[cidx] * (1. - f) + gastodust[didx] * f;
        f = (rad - rvals[aidx]) / (rvals[cidx] - rvals[aidx]);
        *gtd = A * (1. - f) + B * f;
        
        // Check to make sure this isn't below the default value.
        // This is a little trickier as you want it to vary. Think 
        // about what value you want here.
        if (*gtd < 50.){
	        *gtd = 50.;
        }

    }

}

/******************************************************************************/

void
velocity(double x, double y, double z, double *velocity){

    // Assume standard Keplerian rotation. 
    velocity[0] = sqrt(6.67e-11 * MSTAR * 2e30 / sqrt(x*x + y*y + z*z))*sin(atan2(y,x));
    velocity[1] = -sqrt(6.67e-11 * MSTAR * 2e30 / sqrt(x*x + y*y + z*z))*cos(atan2(y,x));
    velocity[2] = 0.0;    
    
}

/******************************************************************************/

