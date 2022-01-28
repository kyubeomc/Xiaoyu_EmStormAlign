clc;
close all;
clear ;
tic;
Deform_robust('2-Cep164-647-cntrl_0001_list_cluster1_20150414_30_s');
[sumim,regi,allIm] = SuperResAlignL('2-Cep164-647-cntrl_0001_list_cluster1_20150414_30_s_tran_redu_robust',5,50, 175, 0.41,100,90,3,10); 
toc;