function   Deform_robust(file)
% by Xiaoyu Shi @ Bo Huang Lab, UCSF
%circlulate + remove >1.5 stand deviations
%% input
%   file: the txt file name of the molecule list.     e.g. 'Cep164-647-high_0005_list_auto_cluster_robust_20140926-20150414'
%   see demo for the structure of the molecule list
%% output
%   A txt file of deformed structures

%% load data

clc
addpath(genpath(pwd));
inname = [file,'.txt'];
rawlist = importdata(inname);
rawdata = rawlist.data;
size0 = size(rawdata);
n0 = size0(1);    % n0 = number of molecules in the raw data
nframe = rawdata(n0,13);

%% transform ellipse to circle

a=[];
b=[];
phi=[];
x0=[];
y0=[];
x=[];
y=[];
trredata=[]; %transformed and reduced data

%------use average diameter of ref channel for deformation-------

Rave = AveRadAB(rawdata,nframe);

%---------------------------------------------------------------

for i=1:nframe
   i
   newdata1 = rawdata(rawdata(:,13)==i,:);
   
   x1 =rawdata(rawdata(:,13)==i,4);
   y1 =rawdata(rawdata(:,13)==i,5);
   
% % % % % % % % % % % % % % % % % % %    
    [xc0,yc0,Re]=circfit(x1,y1);
    xc1=x1-xc0;
    yc1=y1-yc0;
    d1=[num2str(Re),'   nm.'];
    fprintf('The average radius of all data is :  %s\n',d1)

    %convert Cartesian coordinates to polar coordinates

    [theta,r]=cart2pol(xc1,yc1);
    Polar=[theta,r];

    % remove outliers  
    I=(r - Re) > 1.5*std(r);
    [newcartesion1,PS]=removerows(newdata1,'ind',I);
    xc02 = newcartesion1(:,4); 
    yc02 = newcartesion1(:,5);
    
    [theta2,r2]=cart2pol(xc02,yc02);
    
    I=r2<0.5*Re;
    [newcartesion,PS]=removerows(newcartesion1,'ind',I);
    xc2 = newcartesion(:,4); 
    yc2 = newcartesion(:,5);
% % % % % % % % % % % % % % % % % % % % %      
   XY = [xc2 yc2];
   A = EllipseFitByTaubin(XY); 
   [G] = AtoG(A);
   %tran = fit_ellipse(xc2,yc2)
   
   phi=[phi G(5)];
   a=[a G(3)];
   b=[b G(4)];
   x0=[x0 G(1)];
   y0=[y0 G(2)];

   %move data
    [xm,ym]=move(x0(i),y0(i),xc2,yc2);

    %rotate data
    [xr,yr]=rotate(-phi(i),xm,ym);

    %transiform ellipse to circle
    [xcir,ycir]=circle(a(i),b(i),xr,yr,Rave);
   
   x = [x;xcir];
   y = [y;ycir];   
   
   trredata = [trredata;newcartesion];
   
end

trredata(:,4)=x;
trredata(:,5)=y;
nr=size(x);
nr0=nr(1);

%% -----------------------write----------------------
filename = [file,'_tran_redu_robust','.txt'];   %transformed to circles
f = fopen(filename,'wt');
  
%Write the header:

    cas = ['Cas',num2str(nr0)];
    header = {cas 'X' 'Y' 'Xc' 'Yc' 'Height' 'Area' 'Width' 'Phi' 'Ax' 'BG' 'I' 'Frame' 'Length' 'Link' 'Valid' 'Z' 'Zc'};
    fprintf(f,'%s\t',header{1:end-1});
    fprintf(f,'%s\n',header{end});

%Write the data:
    for m = 1:nr0
        fprintf(f,'%g\t',trredata(m,1:end-1));
        fprintf(f,'%g\n',trredata(m,end));
    end

fclose(f);
   


