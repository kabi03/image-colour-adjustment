%% 3.1.1 RGB image to HSI
img = imread('mandrill.png');
img_db = im2double(img);
mex rgb_to_hsi.c
output1 = rgb_to_hsi(img_db);
figure(1);
imshow(img_db); title('original');
figure(2);
imshow(output1); title('output');

%% 3.1.2 HSI image to RGB
img = imread('hsi.png');
img_db = im2double(img);
mex hsi_to_rgb.c
output2 = hsi_to_rgb(output1);
figure(1);
imshow(img_db); title('original');
figure(2);
imshow(output2); title('output');
%% 3.1.3 RGB image to YCbCr
img = imread('mandrill.png');
img_db = im2double(img);
mex rgb_to_ycbcr.c
output3 = rgb_to_ycbcr(img_db);
figure(1);
imshow(img_db); title('original');
figure(2);
imshow(output3); title('output');

%% 3.1.4 YCbCr image to RGB
img = imread('ycbcr.png');
img_db = im2double(img);
mex ycbcr_to_rgb.c
output4 = ycbcr_to_rgb(output3);
figure(1);
imshow(img_db); title('original');
figure(2);
imshow(output4); title('output');
%% 3.2.1 Change Hue
img = imread('mandrill.png');
img_db = im2double(img);
mex change_hue.c
output5 = change_hue(img_db,90);
mex hsi_to_rgb.c
output5 = hsi_to_rgb(output5);
figure(1);
imshow(img_db); title('original');
figure(2);
imshow(output5); title('output');
%% 3.2.2 Change Saturation
img = imread('mandrill.png');
img_db = im2double(img);
mex change_saturation.c
output6 = change_saturation(img_db,30);
mex hsi_to_rgb.c
output6 = hsi_to_rgb(output6);
figure(1);
imshow(img_db); title('original');
figure(2);
imshow(output6); title('output');

%% 3.2.3 Apply point transform
clear all;
close all;
clc;
img = imread('mandrill.png');
img_db = im2double(img);
mex apply_point_tfrm.c
output1 = apply_point_tfrm(img_db,1.5,0.1);
figure(1);
imshow(img_db); title('original');
figure(2);
imshow(output1); title('output');

%% 3.2.4 Spatial Filter
clear all;
close all;
clc;
img = imread('mandrill.png');
img_db = im2double(img);
mex spatial_filter.c
H= [1 1 1 ; 0 0 0; -1 -1 -1];
output1 = spatial_filter_colour(img_db,H);
figure(1);
imshow(img_db); title('original');
figure(2);
imshow(output1); title('apply spatial filter');
%% 3.2.5 Histogram Equalize colour
clear all;
close all;
clc;
img = imread('mandrill.png');
mex histogram_equalize.c
output1 = histogram_equalize(img);
figure(1);
imshow(img); title('original');
figure(2);
imshow(output1); title('Histogram Equalize colour');
%% 3.2.6 Histogram Equalize luma
img = imread('mandrill.png');
img_db = im2double(img);
mex histogram_equalize_luma.c
output1 = histogram_equalize_luma(img_db);
figure(1);
imshow(img_db); title('original');
figure(2);
imshow(output1); title('output');