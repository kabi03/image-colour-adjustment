#include "mex.h"
#include <inttypes.h>
#include <math.h>

void hsi_to_rgb(int width, int height, int channels, const double *input,double *output){
    
    int c1 = 0;
    int c2 = width*height;
    int c3 = 2*width*height;
    

    const double *inputH = &input[c1];
    const double *inputS = &input[c2];
    const double *inputI = &input[c3];

    double *out_red = &output[c1];
    double *out_green = &output[c2];
    double *out_blue = &output[c3];
 
    double PI = 3.14159265359;
    double PI2 = 3.14159265359*2;
    double PI2by3 = PI*2/3;
    double PI4by3 = PI*4/3;
    
    for (int c = 0; c <channels; c++){
        
        for(int i = 0 ; i<height; i++){
            for(int j = 0; j<width; j++){

                int ind = j*height + i;
                
                double H = inputH[ind];
                double S = inputS[ind];
                double I = inputI[ind];
                
                double HP = H - PI*2/3;
                double HPP = H - PI4by3; 
                
                if(H >= 0 && H < PI2by3){

                    out_blue[ind] = I*(1-S);
                    out_red[ind] = I*(1+ (S*cos(H))/(cos(PI/3 - H)));
                    out_green[ind] = 3*I - (out_red[ind] + out_blue[ind]);
                    
                } else if(H >= PI2by3 && H < PI4by3){
                    out_red[ind] = I*(1-S);
                    out_green[ind] = I*(1+ S*cos(H-2*PI/3)/cos(PI - H));
                    out_blue[ind] = 3*I - (out_red[ind] + out_green[ind]);
                    
                } else if(H >= PI4by3 && H <= PI2){
                    out_green[ind] = I*(1-S);
                    out_blue[ind] = I*(1+ S*cos(H-4*PI/3)/cos(5*PI/3 - H));
                    out_red[ind] = 3*I - (out_green[ind] + out_blue[ind]);
                }
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    if (nrhs != 1)
        mexErrMsgTxt("Too many outputs");
    
    if (nlhs != 1)
        mexErrMsgTxt("Too many inputs");
    
    const mxArray *img = prhs[0];
    double *img_ptr = (double *)mxGetData(img);
    
    if ( !mxIsDouble(img))
        mexErrMsgTxt("image is not double");
    
    mwSize ndims = mxGetNumberOfDimensions(img);
    const mwSize *dims = mxGetDimensions(img);
    
    int height_img   = dims[0];
    int width_img    = dims[1];

    int channels = (ndims == 2) ? 1 : dims[2];
    mxClassID input_type = mxGetClassID(img);
    
    mxArray *output = mxCreateNumericArray(ndims, dims, input_type, mxREAL);
    plhs[0] = output;
    
    if (mxIsDouble(prhs[0])){
        double *out_ptr = (double *)mxGetData(output);
        hsi_to_rgb(width_img, height_img, channels, img_ptr, out_ptr);    
    }
    
}