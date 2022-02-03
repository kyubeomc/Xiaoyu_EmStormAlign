%% Fit centroid of EM image 
% 2022-02-01
% Kyubeom Choi

clc;
clear;
close all;
addpath(genpath(pwd));

% the input inside dir is the location of the folder holding the data image
% add '\*.fileformat' at the end
% we could make this process like interface
input = '/Users/kyubeom/Desktop/QBUMCODE/Xiaoyu_EmStormAlign/screenshot/Emimage4New.png';


%% finding the centroid and making the coordinate data in to txt file
 % cell type to char vector form
    imdata = imread(input); % reading out image file
    Gimdata = 256 - im2gray(imdata); % make the image value gray scale 
    
    
    Bimdata = imbinarize(Gimdata,0.93); % binarizing image data
    
    gGimdata = medfilt2(Bimdata);
    % finding the centroid
  stat = regionprops(gGimdata,'centroid');



 figure(1); imshow(gGimdata); hold on; 

dataXY = []; 
% plotting the centroid on top of the image to check it works. 
for stati=1: numel(stat)
    plot(stat(stati).Centroid(1),stat(stati).Centroid(2),'ro');
 
end
