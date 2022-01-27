function [aveR,stdR,Rmin,Rmax]=AveRadAB(rawdata_ref0,n_frame) % n_frame is the number of frames.

% output in pixel

clc
addpath C:\Users\linde\Desktop\matlabp % changed this 

R=[];

for i=1:1:n_frame 
    i 
   
	xc =rawdata_ref0(rawdata_ref0(:,13)==i,4); % select frameN i xc
	yc =rawdata_ref0(rawdata_ref0(:,13)==i,5);
     
%     % Get optimized radius
%      [x0,y0,Re]=circfit(xc,yc);

% Get optimized radius: fit to ellips, radius = perim/2/pi;
    XY = [xc yc];
    A = EllipseFitByTaubin(XY); % input : xc,yc in frameN i  output : 
    [G] = AtoG(A); % phi=G(5);a=G(3);b=G(4);x0=G(1);y0=G(2);       G = [Center(1:2), Axes(1:2), Angle]'
    perim = ellipsePerimeter(G'); % finding the peremiter of the ellipse 둘레
    
  
    
    Re=perim/2/pi; % perimeter / (2*pi) 
  % not sure this means.
  
  
    d1=[num2str(Re),'   nm.']; 
    i1=num2str(i)
    fprintf('The radius of fame %s\n is :  %s\n\n',i1,d1) 
    
    R = [R Re];

end
% R = [] includes radias for i frame

[Rmin,Imin] = min(R); % don't actually need Imin
[Rmax,Imax] = max(R);


aveR = mean(R);
stdR = std(R);

% %-------------write---------------
% 
% filename = ['D:\Xiaoyu\MatlabAnalysis\3dalignment\Diameter\DiameterAB_', file,'.txt'];  
% fileID = fopen(filename,'wt');
% fprintf(fileID,'The average diameter of %d frames is :  %d +- %d nm.\n',n_frame,round(aveD),round(stdD));
% fprintf(fileID,'The minimum diameter is frame NO. %d :  %d nm.\n',Imin,round(Dmin));
% fprintf(fileID,'The maxamum diameter is frame NO. %d :  %d nm.\n\n',Imax,round(Dmax));
% fprintf(fileID,'The pixel size is:  %d nm.\n\n',p);
% fprintf(fileID,'diameter\n');
% fprintf(fileID,'%f\n',D);
% 
% 
% fclose(fileID);
% 
% fprintf('The average diameter of %d frames is :  %d +- %d nm.\n',n_frame,round(aveD),round(stdD));
% fprintf('The minimum diameter is frame NO. %d :  %d nm.\n',Imin,round(Dmin));
% fprintf('The minimum diameter is frame NO. %d :  %d nm.\n',Imax,round(Dmax));