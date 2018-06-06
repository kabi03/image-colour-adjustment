#include "mex.h"
#include <inttypes.h>
#include <math.h>


void rgb_to_ycbcr(int wid, int heig, int chls, const double *input,double *output) {
    
    int c1 = 0;
    int c2 = wid*heig;
    int c3 = 2*wid*heig;
        

    const double *input_red = &input[c1];
    const double *input_green = &input[c2];
    const double *input_blue = &input[c3];

    double *Y = &output[c1];
    double *Cb = &output[c2];
    double *Cr = &output[c3];
    
    for (int c = 0; c <chls; c++){ 
        for(int i = 0 ; i<heig; i++){
            for(int j = 0; j<wid; j++){
                int ind = j*heig + i;
                
                double R = input_red[ind];
                double G = input_green[ind];
                double B = input_blue[ind];
                
                Y[ind]  = 0.299*R + 0.587*G + 0.144*B;
                Cb[ind] = -0.168736*R  - 0.331264*G + 0.5*B + 0.5;
                Cr[ind] = 0.5*R - 0.418688*G - 0.081312*B + 0.5;   
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
        rgb_to_ycbcr(width_img, height_img, channels, img_ptr, out_ptr);    
    }
    
}
