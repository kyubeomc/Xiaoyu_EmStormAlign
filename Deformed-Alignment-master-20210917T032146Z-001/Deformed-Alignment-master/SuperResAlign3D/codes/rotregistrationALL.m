

function [aligniter,iterIm,s,n_frame,allIm,data0] = rotregistrationALL(images,data0,usfac,angrange,angstep,n_iteration,zoomfactor)
% by Xiaoyu Shi @ Bo Huang Lab, UCSF

imagesize = size(images);
n_frame = imagesize(3);
s = imagesize(1);
sumim = zeros(s,s); 
phase = [];
CCall = [];
aligniter = [];

%% superimpose all images
for index = 1:n_frame
    sumim = sumim + images(:,:,index);
end

%% align images for n_iterations iterations
for iter = 1:n_iteration
    
        fprintf('--------iteration %d --------\n',iter);
        
        % align image k, the reference is all the other images
        for k=1:n_frame;
            
            fprintf('image %d\n',k);
            im2 = images(:,:,k);
            im1 = sumim-images(:,:,k);
            ft1 = fft2(im1);
            ft2 = fft2(im2);

            [Greg,row_shift,col_shift,Angle] = rotregistration(ft1,ft2,usfac,angrange,angstep);
                    
            alignk = [Angle row_shift col_shift];   %angle row_shift col_shift 
       
            A = abs(ifft2(Greg));
            images(:,:,k)= A;
            allIm(:,:,k,iter)=A;
            aligniter(k,:,iter)= alignk;

        end
       
        [sumim]= sumimage(s,n_frame,images);
        iterIm(:,:,iter)=sumim;
        
%% change molecule list
    x2 = [];
    y2 = [];

    for i=1:n_frame
        newdata0 = data0(data0(:,13)==i,:);
        x0 = newdata0(:,4);
        y0 = newdata0(:,5);
        
        da0 = aligniter(i,1,iter)/180*pi;
        dx0 = aligniter(i,3,iter)/zoomfactor;
        dy0 = aligniter(i,2,iter)/zoomfactor;

        %transform                    

        x1 = x0*cos(-da0) - y0*sin(-da0)+ dx0; 
        y1 = x0*sin(-da0) + y0*cos(-da0)+ dy0;  

        x2 = [x2;x1];
        y2 = [y2;y1];

    end

    data0(:,4)=x2;
    data0(:,5)=y2;
      
end
 
return

    
