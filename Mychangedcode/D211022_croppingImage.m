%% cropping Image using function readPoints
% Kyubeom Choi
% 2021-10-14

clear all;
clc;
close all;
% Read data
FileName = 'Cep30_30';
Format = '.png';
Name = [FileName,Format];
imdata = imread(Name);


% make the image gray image
imdataG = rgb2gray(imdata);
figure(1); imshow(imdataG);

% select 2 points in the imdata to collect coordinates
ImageCor = readPoints(imdataG);

% cropped image's width and heidth
ImageWid = ImageCor(1,2) - ImageCor(1,1);
ImageHeid = ImageCor(2,2) - ImageCor(2,1);
ImageCor_Dif = [ImageCor(:,1)',ImageWid, ImageHeid];


CroppedImage = imcrop(imdataG,ImageCor_Dif);


[CroppedImageSize_x,CroppedImageSize_y] = size(CroppedImage); 

% if the cropped image coordinate is a odd number add 1 to it
    if mod(400-CroppedImageSize_x,2) == 1
        CroppedImageSize_x = CroppedImageSize_x + 1;
    end
    
    if mod(400-CroppedImageSize_y,2) == 1
        CroppedImageSize_y = CroppedImageSize_y + 1;
    end

x = (400-CroppedImageSize_x)/2;
y = (400-CroppedImageSize_y)/2;
PadImage = padarray(CroppedImage,[x,y],0,'both');
figure(3); imshow(PadImage);
NewImageName = [FileName, 'New', Format];
imwrite(PadImage, NewImageName);



function pts = readPoints(image, n)
%readPoints   Read manually-defined points from image
%   POINTS = READPOINTS(IMAGE) displays the image in the current figure,
%   then records the position of each click of button 1 of the mouse in
%   figure, and stops when another button is clicked. The track of points
%   is drawn as it goes along. The result is a 2 x NPOINTS matrix; each
%   column is [X; Y] for one point.
% 
%   POINTS = READPOINTS(IMAGE, N) reads up to N points only.
    if nargin < 2
        n = Inf;
        pts = zeros(2, 0);
    else
        pts = zeros(2, n);
    end
imshow(image); title('select 2 points diagnally');   % display image
xold = 0;
yold = 0;
k = 0;
hold on;           % and keep it there while we plot
    while 1
        [xi, yi, but] = ginput(1);      % get a point
        if ~isequal(but, 1)             % stop if not button 1
            break
        end
        k = k + 1;
        pts(1,k) = xi;
        pts(2,k) = yi;
          if xold
              plot([xold xi], [yold yi], 'go-');  % draw as we go
          else
              plot(xi, yi, 'go');         % first point on its own
          end
          
          if isequal(k, n)
              break
          end
          xold = xi;
          yold = yi;
      end
hold off;

    if k < size(pts,2)
        pts = pts(:, 1:k);
    end
end
