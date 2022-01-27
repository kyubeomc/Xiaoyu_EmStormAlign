function [sumim,sumimregi,allIm]= SuperResAlign3D2C(file,imsize,zoomfactor, pixelsize,...
    photonpercount,angrange,angstep,rotrange,rotstep,n_iteration)
%% by Xiaoyu Shi @ Bo Huang Lab, UCSF   
%   2C 3D alignment
%% input
%   file is the name of the txt file of molecule list.     e.g. 'Cep164-647-high_0005_list_auto_cluster_robust_20140926-20150414'
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
%   The output is a txt file of molecule list of aligned clusters.  

%   Xiaoyu Shi@Bo Huang Lab at UCSF
%   Apr 11, 2018

%% load coordinates
tic; 
addpath(genpath(pwd));
fprintf('Loading your molecule list ...\n');
sortdata(file);
filename = [file,'_s','.txt'];

[MList1,data0] = LoadMTxtList2C(filename);

fprintf('Generating image array ...\n');
HRdimension = 3; % 3D array
[imagesRaw] = NormalizedGaussianHR( MList1,imsize,zoomfactor, pixelsize, photonpercount, HRdimension );
nframe = size(imagesRaw,3);
s=size(imagesRaw,1);
images0 = imagesRaw;
for i=1:1:nframe
    [images0(:,:,i)]= LowerImages(imagesRaw(:,:,i),10);   % desity weight, lower contrast
end

[sumim]= sumimage(s,nframe,images0);    

%% align & save aligned molec list
fprintf('Aligning your super-res images ...\n');
[aligniter,sumimregi,s,allIm,data,data5] = rotregistrationALL3D2C(images0,data0,angrange,angstep,n_iteration,rotrange,rotstep);


%% write
writelist11(data5,file);  % output an aligned molecule list

%% display

figure;
subplot(1,2,1);imshow(autocontrast(sumim));title('to be registered')
subplot(1,2,2);imshow(autocontrast(sumimregi));title('registered')

toc;

