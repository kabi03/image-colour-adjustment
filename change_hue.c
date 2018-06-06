#include "mex.h"
#include <inttypes.h>
#include <math.h>

double findMin(double x, double y, double z){
    double min = 1;
    if (x<1)
        min = x;
    if (y<min)
        min = y;
    if (z<min)
        min = z;
    return min;
}

void change_hue(int wid, int heig, int chls, const double *input,double *output, const int hue){
    
    int c1 = 0;
    int c2 = wid*heig;
    int c3 = 2*wid*heig;
    
    
    const double *input_red = &input[c1];
    const double *input_green = &input[c2];
    const double *input_blue = &input[c3];
    
  
    double *H = &output[c1];
    double *S = &output[c2];
    double *I = &output[c3];
  
    for (int c = 0; c <chls; c++){
        for(int i = 0 ; i<heig; i++){
            for(int j = 0; j<wid; j++){
                int ind = j*heig + i;
                
                double R = input_red[ind];
                double G = input_green[ind];
                double B = input_blue[ind];
                
                double numer = 0.5*((R - G) + (R - B));
                double denom = sqrt(pow((R - G),2)+ (R - B)*(G - B));
                double theta = acos(numer/denom);
                double cmin;
                
                if(B <= G){
                    H[ind] = theta;
                }
                if(B > G){
                    H[ind] =2*(3.141592)-theta;
                }

                int H_deg = H[ind]*180/(3.141592);
                H_deg = H_deg + hue;
                H_deg = H_deg%360;
                H_deg = (H_deg + 360) % 360;
                H[ind] = H_deg*(3.141592)/180;
                
                
                cmin = findMin(R,G,B);

                S[ind] = 1 - cmin*(3/(R+G+B));
               
                I[ind] = (R+G+B)/3;
            }
        }    
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    if (nrhs != 2)
        mexErrMsgTxt("Too manny");
    
    if (nlhs != 1)
        mexErrMsgTxt("Too many outputs");
    
    const mxArray *img = prhs[0];
    int hue = mxGetScalar(prhs[1]);
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
    
    if (mxIsDouble(prhs[0]))
    {
        double *out_ptr = (double *)mxGetData(output);
        change_hue(width_img, height_img, channels, img_ptr, out_ptr, hue);    
    }    
}