
function [aligniter,sumim,s,allIm,data3] = rotregistrationALL3D(images,data0,angrange,angstep,n_iteration,rotrange,rotstep,imsize,zoomfactor, pixelsize, photonpercount)

imagesize = size(images);
n_frame = imagesize(3);
s = imagesize(1);
sumim = zeros(s,s); 
phase = [];
CCall = [];
aligniter = [];

%% superimpose all images
[sumim]= sumimage(s,n_frame,images);  

        %figure; imshow(autocontrast(sumim));
%% split master molecule list to frames
alldataf = struct([]);
for i=1:n_frame
    dataff = struct('data',data0(data0(:,13)==i,:));
    alldataf = [alldataf dataff];
end

%% align images for n_iterations iterations
iterIm = [];
for iter = 1:n_iteration
    
        fprintf('--------iteration %d --------\n',iter);
     
        % align image k, the reference is all the other images
        for k=1:n_frame;
            %k
            fprintf('image %d\n',k);
            im2 = images(:,:,k);
            im1 = sumim-images(:,:,k);
            ft1 = fft2(im1);
            dataf0 = alldataf(1,k).data;

            [allim2x,dataf3,alignk] = rotregistration3D(ft1,dataf0,angrange,angstep,rotrange,rotstep,imsize,zoomfactor, pixelsize, photonpercount);
            %[Mlist3,Anglex,Anglez,row_shiftx,col_shiftx] = rotregistration3D(ft1,ft2,usfac,angrange,angstep,Mlist2,rotrange,rotstep);         
            % alignk = [Anglex Anglez row_shift col_shift];  
       
            A = abs(allim2x);
            images(:,:,k)= A;
            allIm(:,:,k,iter)=A;
            aligniter(k,:,iter)= alignk;      
            alldataf(1,k).data = dataf3;

        end
       
        [sumim]= sumimage(s,n_frame,images);
        iterIm(:,:,iter)=sumim;
        
        % change molecule list
        data3 = [];
        for m = 1:n_frame
            data33 = alldataf(1,m).data;
            data3 = [data3;data33];
        end 

%       
end
 
return

    
