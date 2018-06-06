#include <stdio.h>
#include <mex.h>
#include <math.h>

void spatial_filter_colour(double *img , double *output, int img_height, int img_width, int img_channels, double contrast, double brightness){
    double input3D[1000][1000][3] = {0};
    double r[1000][1000] = {0};
    double g[1000][1000] = {0};
    double b[1000][1000] = {0};
    
 
    for(int x=0; x<img_height; x++){
        for(int y=0; y<img_width; y++){
            for(int z=0; z<img_channels; z++){
                int ind = x+y*img_height+z*(img_height*img_width);
                
                if (z ==0)
                    r[x][y] =img[ind];
                else if (z==1)
                    g[x][y] =img[ind];
                else if (z==2)
                    b[x][y] =img[ind];
            }
        }
    }
    
    for (int x = 0; x < img_height; x++){
        for (int y = 0; y < img_width; y++){
            r[x][y] =contrast*r[x][y] +brightness;
            b[x][y] =contrast*b[x][y] +brightness;
            g[x][y] =contrast*g[x][y] +brightness;
            
            if (r[x][y] >= 1)
                r[x][y] = 1;
            else if (r[x][y] <= 0)
                r[x][y] = 0;
            else
                continue;
            
            if (b[x][y] >= 1)
                b[x][y] = 1;
            else if (b[x][y] <= 0)
                b[x][y] = 0;
            else
                continue;
            
            if (g[x][y] >= 1)
                g[x][y] = 1;
            else if (g[x][y] <= 0)
                g[x][y] = 0;
            else
                continue;
        }
    }
    
    for(int x=0; x<img_height; x++){
        for(int y=0; y<img_width; y++){
            for(int z=0; z<img_channels; z++){
                int ind = x+y*img_height+z*(img_height*img_width);
                if (z ==0)
                    output[ind]= r[x][y];
                else if (z==1)
                    output[ind]= g[x][y];
                else if (z==2)
                    output[ind]= b[x][y];
            }
        }
    }
    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    double contrast, brightness;
    
    mxArray *img_m;
    mxArray *output_m;
    
    double *img;
    double *output;
    
    
    img_m = mxDuplicateArray(prhs[0]);
    

    const mwSize *image_dims = mxGetDimensions(prhs[0]);
    int img_height = image_dims[0];
    int img_width = image_dims[1];
    int img_channels = image_dims[2];
    
    contrast = mxGetScalar(prhs[1]);
    brightness = mxGetScalar(prhs[2]);
    
   
    output_m = plhs[0] = mxCreateNumericArray (3, image_dims, mxDOUBLE_CLASS, mxREAL );
    

    img = (double*)mxGetPr(img_m);
    output = (double*)mxGetPr(output_m);
    
    spatial_filter_colour(img, output, img_height, img_width, img_channels, contrast, brightness);
    
    return;
}