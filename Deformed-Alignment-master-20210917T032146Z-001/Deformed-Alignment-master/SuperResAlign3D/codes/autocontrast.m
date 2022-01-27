function [im1]= autocontrast(im0)
% by Xiaoyu Shi @ Bo Huang Lab, UCSF
inputImageMinVal = min(im0(:));
inputImageMaxVal = max(im0(:));

im1 = (im0 - inputImageMinVal) / (inputImageMaxVal - inputImageMinVal);