function [sumim]= sumimage(s,n_frame,images)
% by Xiaoyu Shi @ Bo Huang Lab, UCSF
sumim = zeros(s,s);
for index = 1:n_frame
    sumim = sumim + images(:,:,index);
end