%% finding the centroid of the image and making as our center
% 2021-10-16 : find the centroid in side the 
% Kyubeom Choi
% Image segmentation using threshold and histogram

clc
clear all;
close all;
addpath(genpath(pwd));
imdata = imread('/Users/kyubeom/Desktop/QBUMCODE/Xiaoyu_EmStormAlign/screenshot/Emimage1.png'); % reading out image file
Gimdata = im2gray(imdata); % make the image value gray scale! 

Bimdata = imbinarize(Gimdata);

%% making the image binary and find the centroid of the images. 
% This way we can find xc and yc more easy 




% we could use histogram to check the threshold ? and make the code
% adjustable

stat = regionprops(Bimdata,'centroid'); % finding the centroid of the BWimage this makes easier to find the center of the picture
figure(1); imshow(imdata); hold on; 
dataXY = []; 
% plotting the centroid on top of the image to check it works. 
for stati=1: numel(stat)
    plot(stat(stati).Centroid(1),stat(stati).Centroid(2),'r.');
    dataXY = [dataXY;stat(stati).Centroid(1),stat(stati).Centroid(2)];
end


writematrix(dataXY,'tabledata.txt'); % saving the martix in txt file. 

matrixdata = importdata('tabledata.txt'); % reading txt to array. 





