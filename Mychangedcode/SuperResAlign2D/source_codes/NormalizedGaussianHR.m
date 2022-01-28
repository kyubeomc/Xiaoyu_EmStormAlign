function [ HRimage ] = NormalizedGaussianHR( MList,imsize,zoomfactor, pixelsize, photonpercount, HRdimension )
%   NormalizedGaussianHR function by Yina Wang
%   Rendering Gaussian width: molecule width / sqrt(photon number)
%   input:
%   Mlist: Insight3 molecule list
%   imsize:  x and y image size, assuming squred image
%   zoomfactor: zoomed factor of the HR image relative to raw data
%   pixelsize: in nm
%   photonpercount: photon number per grey count
%   HRdimension: output 2D or 3D data array. 

coords(:,1) = MList.Xc+imsize/2;
coords(:,2) = MList.Yc+imsize/2;
coords(:,3) = MList.I;
coords(:,4) = MList.Width;
coords(:,5) = MList.Frame;

% Offset and rescale positions
coords = single(coords);
coords(:,1:2) = coords(:,1:2)*zoomfactor;
coords(:,1:2) = coords(:,1:2) + 0.5; % offset to the pixel center

% Filter localizations outside the FOV
binimsize = ceil(imsize*zoomfactor);
keep = coords(:,1)>=0;
keep = keep & coords(:,2)>=0;
keep = keep & coords(:,1)<=binimsize;
keep = keep & coords(:,2)<=binimsize;
coords = coords(keep,:);
nframe = max(coords(:,5));

% grid for gaussian rendering
gridsize = round(max(coords(:,4))/pixelsize * zoomfactor);
[X,Y] = meshgrid(-gridsize:gridsize,-gridsize:gridsize);
gaussiankernel = @(x,y,w) exp(-((X-x).^2+(Y-y).^2)/(2*w^2))/sqrt(2*pi*w^2);

% rendering localizations
HRimage = [];
if HRdimension == 2
    HRimage = zeros(binimsize+2*gridsize,binimsize+2*gridsize);
    for i = 1:size(coords,1)
        x0 = coords(i,1) + gridsize;
        y0 = coords(i,2) + gridsize;
        w = (coords(i,4) / sqrt(coords(i,3)*photonpercount)) / pixelsize * zoomfactor;
        
        x = x0 - (floor(x0) + 0.5); 
        y = y0 - (floor(y0) + 0.5); 
        molim = gaussiankernel(x,y,w);
        
        HRimage(ceil(y0)-gridsize:ceil(y0)+gridsize,ceil(x0)-gridsize:ceil(x0)+gridsize) = ...
            HRimage(ceil(y0)-gridsize:ceil(y0)+gridsize,ceil(x0)-gridsize:ceil(x0)+gridsize) + molim;
    end
    HRimage = HRimage(gridsize+1:end-gridsize,gridsize+1:end-gridsize);
%        %lower contrast
%     [HRimage] = LowerImages(HRimage,3);
elseif HRdimension == 3
    HRimage = zeros(binimsize+2*gridsize,binimsize+2*gridsize,nframe);  
    for i = 1:size(coords,1)
        x0 = coords(i,1) + gridsize;
        y0 = coords(i,2) + gridsize;
        w = (coords(i,4) / sqrt(coords(i,3)*photonpercount)) / pixelsize * zoomfactor;
        
        x = x0 - (floor(x0) + 0.5); 
        y = y0 - (floor(y0) + 0.5); 
        molim = gaussiankernel(x,y,w);
        
        HRimage(ceil(y0)-gridsize:ceil(y0)+gridsize,ceil(x0)-gridsize:ceil(x0)+gridsize, coords(i,5)) = ...
            HRimage(ceil(y0)-gridsize:ceil(y0)+gridsize,ceil(x0)-gridsize:ceil(x0)+gridsize, coords(i,5)) + molim;
%         %lower contrast
%         Image0=HRimage(ceil(y0)-gridsize:ceil(y0)+gridsize,ceil(x0)-gridsize:ceil(x0)+gridsize, coords(i,5)); 
%         [HRimage(ceil(y0)-gridsize:ceil(y0)+gridsize,ceil(x0)-gridsize:ceil(x0)+gridsize, coords(i,5)) ] = LowerImages(Image0,3);
    end
    HRimage = HRimage(gridsize+1:end-gridsize,gridsize+1:end-gridsize,:);
else
    HRimage = [];
    
end

