function [sumim,regi,allIm]= SuperResAlignL(file,imsize,zoomfactor, pixelsize, photonpercount,usfac,angrange,angstep,n_iteration)
% by Xiaoyu Shi @ Bo Huang Lab, UCSF

%% input
%   file is the name of the txt file of molecule list.     e.g. 'Cep164-647-high_0005_list_auto_cluster_robust_20140926-20150414'
%           see demo for the structure of the molecule list
%   imsize:  x and y image size, assuming squred image, in pixel   !! make
%                           sure the imsize is bigger than your structure!!
%   zoomfactor: zoomed factor of the HR image relative to raw data
%   pixelsize: in nm
%   photonpercount: photon number per grey count
%   usfac is the up ssampling for alignment. Higher number does not extends the computing time.
%   angrange is the range you wanted to scan your angle for alignment
%   angstep is the interval angle
%   n_iteration is the total number of the iterations for alignment.

%   You should start to see good result from iteration 3.  The cross correlation converges
%   after 10 iterations usually.

%% output
%   A txt file of molecule list of aligned clusters.  
%   sumim: average of original iamges.
%   regi:  average of aligned iamges.
%   allIm: alignment information.

%   Xiaoyu Shi@Bo Huang Lab at UCSF
%   Apr 11, 2018

%% load coordinates
tic; 
addpath(genpath(pwd));

filename = [file,'.txt'];
fprintf('Loading your molecule list ...\n');
[MList,data0]=LoadMTxtList(filename);

fprintf('Generating image array ...\n');
HRdimension = 3; % 3D array
[imagesRaw] = NormalizedGaussianHR( MList,imsize,zoomfactor, pixelsize, photonpercount, HRdimension );
nframe = size(imagesRaw,3);
images0 = imagesRaw;
for i=1:1:nframe
    [images0(:,:,i)]= LowerImages(imagesRaw(:,:,i),10);   % desity weight, lower contrast
end
    

%% align & save aligned molec list
fprintf('Aligning your super-res images ...\n');
[aligniter,iterIm,s,n_frame,allIm,data] = rotregistrationALL(images0,data0,usfac,angrange,angstep,n_iteration,zoomfactor);

%% write
writelist11(data,file);    % Outputs a txt file of molecule list of aligned clusters.  

%% display
sumim = zeros(s,s);
for index = 1:n_frame
    sumim = sumim + images0(:,:,index);
end
% 가장 마지막 iteration 이 가장 확실한 건가보네???
regi=iterIm(:,:,n_iteration);

AAA = zeros(1,n_iteration)
for index = 1:n_iteration-1
    AAA(index) = sum(sum(((iterIm(:,:,index+1) - iterIm(:,:,index)).^2)))
end



figure;
subplot(1,2,1);imshow(autocontrast(sumim));title('to be registered')
subplot(1,2,2);imshow(autocontrast(regi));title('registered')

toc;

