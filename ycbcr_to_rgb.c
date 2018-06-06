#include "mex.h"
#include <inttypes.h>
#include <math.h>



void ycbcr_to_rgb(int wid, int heig, int chls, const double *input,double *output) {
    
 
    int c1 = 0;
    int c2 = wid*heig;
    int c3 = 2*wid*heig;
    
   
    
    const double *colour_inY = &input[c1];
    const double *colour_inCb = &input[c2];
    const double *colour_inCr = &input[c3];
    

    double *colour_outR = &output[c1];
    double *colour_outG = &output[c2];
    double *colour_outB = &output[c3];
    
    
   
    for (int c = 0; c <chls; c++){  
        for(int i = 0 ; i<heig; i++){
            for(int j = 0; j<wid; j++){
                int ind = j*heig + i;
                
              
                double Y = colour_inY[ind];
                double Cb = colour_inCb[ind];
                double Cr = colour_inCr[ind];
                
            
                colour_outR[ind] = Y +  0*Cb + 1.402*(Cr-0.5);
                colour_outG[ind] = Y -0.344136*(Cb-0.5)  - 0.714136*(Cr-0.5);
                colour_outB[ind] = Y + 1.772*(Cb-0.5) - 0*Cr;
            }   
        } 
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    if (nrhs != 1)
        mexErrMsgTxt("ycbcr_to_rgb can only accept two input argument");
    
    if (nlhs != 1)
        mexErrMsgTxt("ycbcr_to_rgb requires one output argument");
    
    
    
    const mxArray *img = prhs[0];
    
    
    double *img_ptr = (double *)mxGetData(img);
    
    
    if ( !mxIsDouble(img))
        mexErrMsgTxt("ycbcr_to_rgb can only accept images of type 'double'");


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

        int N = nrhs;

        ycbcr_to_rgb(width_img, height_img, channels, img_ptr, out_ptr);

    }
}