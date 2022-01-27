function [allim2x,dataf3,alignk] = rotregistration3D(buf1ft,dataf0,angrange,angstep,rotrange,rotstep,...
    imsize,zoomfactor, pixelsize, photonpercount)
%% by Xiaoyu Shi @ Bo Huang Lab, UCSF
%%
usfac = 100;

dataf2=dataf0;
dataf3=dataf0;
x0=dataf0(:,4);
y0=dataf0(:,5);
z0=dataf0(:,18)/175;

for l=1:(rotrange/rotstep+1)  % 3d rotation around x
    %l
    thetax=(rotrange/2*(-1)+(l-1)*rotstep)/180*pi;
    x1 = x0;
    y1 = y0*cos(thetax) - z0*sin(thetax);
    z1 = y0*sin(thetax) + z0*cos(thetax);            

    dataf2(:,4) = x1;
    dataf2(:,5) = y1;
    dataf2(:,18) = z1/175;
    xyz(:,:,l) = [x1 y1 z1];

    [ im2 ] = NormalizedGaussian(dataf2,imsize,zoomfactor, pixelsize, photonpercount);
    [im2]= LowerImages(im2,10);
    
            ccc = [];       %%%%%%%%%%%%%%%%%%%%%%%%%%%
            shiftrow = [];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            shiftcol = [];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            diffphase = [];

            for j=1:angrange/angstep  % 3d rotation around z
                %j
                buf2ft = fft2(imrotate(im2,(j-1)*angstep,'crop'));
                % First upsample by a factor of 2 to obtain initial estimate
                % Embed Fourier data in a 2x larger array
                [m,n]=size(buf1ft);
                mlarge=m*2;
                nlarge=n*2;
                CC=zeros(mlarge,nlarge);
                CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
                    fftshift(buf1ft).*conj(fftshift(buf2ft));

                % Compute crosscorrelation and locate the peak 
                CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
                [max1,loc1] = max(CC);
                [max2,loc2] = max(max1);
                rloc=loc1(loc2);cloc=loc2;
                CCmax=CC(rloc,cloc);


                % Obtain shift in original pixel grid from the position of the
                % crosscorrelation peak 
                [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
                        if rloc > md2 
                            row_shift = rloc - m - 1;
                        else
                            row_shift = rloc - 1;
                        end
                        if cloc > nd2
                            col_shift = cloc - n - 1;
                        else
                            col_shift = cloc - 1;
                        end
                        row_shift=row_shift/2;
                        col_shift=col_shift/2;

                        % If upsampling > 2, then refine estimate with matrix multiply DFT

                            %%% DFT computation %%%
                            % Initial shift estimate in upsampled grid
                            row_shift = round(row_shift*usfac)/usfac; 
                            col_shift = round(col_shift*usfac)/usfac;     
                            dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
                            % Matrix multiply DFT around the current shift estimate
                            CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
                                dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);
                            % Locate maximum and map back to original pixel grid 
                            [max1,loc1] = max(CC);   
                            [max2,loc2] = max(max1); 
                            rloc = loc1(loc2); cloc = loc2;
                            CCmax = CC(rloc,cloc);
                            rg00 = dftups(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
                            rf00 = dftups(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);  
                            rloc = rloc - dftshift - 1;
                            cloc = cloc - dftshift - 1;
                            row_shift = row_shift + rloc/usfac;
                            col_shift = col_shift + cloc/usfac;    


                        error = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
                        error = sqrt(abs(error));
                        diffphase=atan2(imag(CCmax),real(CCmax));
                        % If its only one row or column the shift along that dimension has no
                        % effect. We set to zero.
                        if md2 == 1,
                            row_shift = 0;
                        end
                        if nd2 == 1,
                            col_shift = 0;
                        end

                        ccc = [ccc max2];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        shiftrow = [shiftrow row_shift];    %%%%%%%%%%%%%%%%%%%%%
                        shiftcol = [shiftcol col_shift];   %%%%%%%%%%%%%%%%%%%%%
                        shiftphase(:,:,j) = diffphase;     %%%%%%%%%%%%%%%%%%%%%

                        diffphase = [];
                        CC = [];

            end  

            [maxall,Irot] = max(ccc);
            %Irot
            Angle = (Irot-1)*angstep;
            im2 = imrotate(im2,Angle,'crop');
            %buf2ft = fft2(imrotate(im2,Angle,'crop'));
            row_shift = shiftrow(Irot);
            col_shift = shiftcol(Irot);
            diffphase = shiftphase(:,:,Irot);

            allmax(l)=max(ccc);
            allphase(:,:,l)= diffphase;
            allangle(l)= Angle;
            allrowshift(l)= row_shift;
            allcolshift(l)= col_shift;
            
            allim2(:,:,l)=im2;
            

end
[allallmaxl,Ix] = max(allmax);
%Ix
Anglez = allangle(Ix)*(-1);
Anglex = (Ix-1)*rotstep*(-1);
buf2ftx = fft2(allim2(Ix));
row_shiftx = allrowshift(Ix);
col_shiftx = allcolshift(Ix);
diffphasex = allphase(:,:,Ix);
allim2x = allim2(:,:,Ix);

x2 = xyz(:,1,Ix);
y2 = xyz(:,2,Ix);
z2 = xyz(:,3,Ix);

x3 = x2*cos(Anglez/180*pi) - y2*sin(Anglez/180*pi)+ col_shiftx/20; 
y3 = x2*sin(Anglez/180*pi) + y2*cos(Anglez/180*pi)+ row_shiftx/20;

dataf3(:,4) = x3;
dataf3(:,5) = y3;
dataf3(:,18) = z2*175;

alignk = [Anglex Anglez row_shift col_shift];  


%% Compute registered version of buf2ft
    [nr,nc]=size(buf2ftx);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ftx.*exp(i*2*pi*(-row_shiftx*Nr/nr-col_shiftx*Nc/nc));
    Greg = Greg*exp(i*diffphasex);

return

function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the 
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc]=size(in);
% Set defaults
if exist('roff')~=1, roff=0; end
if exist('coff')~=1, coff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
kernr=exp((-i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
return
