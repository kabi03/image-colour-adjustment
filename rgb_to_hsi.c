#include "mex.h"
#include <inttypes.h>
#include <math.h>
#define pi 3.14159265359

double findMin(double a, double b, double c){
    double min = 1;
    if (a<1)
        min = a;
    if (b<min)
        min = b;
    if (c<min)
        min = c;
    return min;
}

void rgb_to_hsi(int wid, int heig, int channels, const double *input,double *output){
    
    int c1 = 0;
    int c2 = wid*heig;
    int c3 = 2*wid*heig;
    
    
    const double *input_red = &input[c1];
    const double *input_green = &input[c2];
    const double *input_blue = &input[c3];
    
    
    double *H = &output[c1];
    double *S = &output[c2];
    double *I = &output[c3];
    
    for (int c = 0; c <channels; c++){
        
        for(int i = 0 ; i<heig; i++){
            for(int j = 0; j<wid; j++){

                int ind = j*heig + i;
                
                double R = input_red[ind];
                double G = input_green[ind];
                double B = input_blue[ind];
                
                double theta = acos((0.5*((R - G) + (R - B)))/(sqrt(pow((R - G),2)+ (R - B)*(G - B))));
                double cmin;
                
                if(B <= G)
                    H[ind] = theta;
                if(B > G)
                    H[ind] =2*pi - theta;
                
                cmin = findMin(R,G,B);
                S[ind] = 1-cmin*(3/(R+G+B));
                I[ind] = (R+G+B)/3; 
            }
        } 
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    const mxArray *img = prhs[0];
    double *img_ptr = (double *)mxGetData(img);
    
    mwSize ndims = mxGetNumberOfDimensions(img);
    const mwSize *dims = mxGetDimensions(img);
    
    int height_img   = dims[0];
    int width_img    = dims[1];
    
    int channels = (ndims == 2) ? 1 : dims[2];
    mxClassID input_type = mxGetClassID(img);
    
    mxArray *output = mxCreateNumericArray(ndims, dims, input_type, mxREAL);
    plhs[0] = output;
    

        double *out_ptr = (double *)mxGetData(output);
        rgb_to_hsi(width_img, height_img, channels, img_ptr, out_ptr);    
    
}