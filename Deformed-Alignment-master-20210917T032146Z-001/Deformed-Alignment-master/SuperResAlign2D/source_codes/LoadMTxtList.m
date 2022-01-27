function [MList,data0] = LoadMTxtList(filename)
% by Xiaoyu Shi @ Bo Huang Lab, UCSF
MListData = importdata(filename);
data0 = MListData.data;
colheaders = genvarname(MListData.colheaders);
colheaders{1} = 'C';
for i = 1:length(colheaders)
    MList.(colheaders{i}) = MListData.data(:, i);
end