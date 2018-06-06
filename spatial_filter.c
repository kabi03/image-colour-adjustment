#include "mex.h"
#include <inttypes.h>
#include <math.h>

void spatial_filter(int wid, int heig, int chls, int wid_msk, int heig_msk,const double *input,const double *mask, double *output) {
    
    for (int c = 0; c < chls; c++){
        
        int colour = c*wid*heig;

        const double *gray = &input[colour];
        const double *msk = &mask[0];
        double *grayout = &output[colour];
        
        int mask_wr = floor(wid_msk/2);
        int mask_hr = floor(heig_msk/2);
        int mask_wr_var = mask_wr;
        int wid_msk_var = 1;
        double temp[wid_msk*heig_msk];
        
        int heig_edge = heig - mask_hr;
        int wid_edge = wid- mask_wr ;
        
        for (int y = mask_hr; y < heig_edge; y++){
            for (int x = mask_wr; x < wid_edge; x++){
                
                int ind = x*heig + y;
                int ind_msk = mask_wr * heig_msk + mask_hr;
                int k = ind_msk - mask_wr * heig_msk;
                for(int j = ind - mask_wr * heig; j <= ind + mask_wr * heig;j = j + heig){
                    
                    if(k <= ind_msk + mask_hr * heig_msk){
                        for (int l = -mask_hr; l <= mask_hr; l++){

                            grayout[ind] = grayout[ind]+ (gray[j+l] * msk[k+l]);
                             
                        }
                        k = k + heig_msk;                        
                    }
                }                
            }                     
        }
    }    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    if (nrhs != 2)
        mexErrMsgTxt("can only accept two input argument");
    
    if (nlhs != 1)
        mexErrMsgTxt("can only accept one output argument");
    
    const mxArray *img = prhs[0];
    const mxArray *mask = prhs[1];
    
    double *img_ptr = (double *)mxGetData(img);
    double *mask_ptr = (double *)mxGetData(mask);
    
    if ( !mxIsDouble(img))
        mexErrMsgTxt("can only accept images of type 'double'");
 
    mwSize ndims = mxGetNumberOfDimensions(img);
    const mwSize *dims = mxGetDimensions(img);
    
    mwSize ndims_msk = mxGetNumberOfDimensions(mask);
    const mwSize *dims1 = mxGetDimensions(mask);
    
    int height_img   = dims[0];
    int width_img    = dims[1];
    
    int height_mask   = dims1[0];
    int width_mask    = dims1[1];
    int channels = (ndims == 2) ? 1 : dims[2];
    
    mxClassID input_type = mxGetClassID(img);
    
    mxArray *output = mxCreateNumericArray(ndims, dims, input_type, mxREAL);
    plhs[0] = output;
    
    if (mxIsDouble(prhs[0])){
        
        double *out_ptr = (double *)mxGetData(output);
        int N = nrhs;
        
        spatial_filter(width_img, height_img, channels, width_mask, height_mask, img_ptr, mask_ptr, out_ptr);
    }
    
}