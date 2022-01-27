function writelist11(data,file)
% by Xiaoyu Shi @ Bo Huang Lab, UCSF

size0 = size(data);
n0 = size0(1);


    %-----------------------write----------------------
    filename = [file,'_aligned','.txt'];   
    f = fopen(filename,'wt');

    %Write the header:

        cas = ['Cas',num2str(n0)];
        header = {cas 'X' 'Y' 'Xc' 'Yc' 'Height' 'Area' 'Width' 'Phi' 'Ax' 'BG' 'I' 'Frame' 'Length' 'Link' 'Valid' 'Z' 'Zc'};
        fprintf(f,'%s\t',header{1:end-1});
        fprintf(f,'%s\n',header{end});

    %Write the data:
        for m = 1:n0
            fprintf(f,'%g\t',data(m,1:end-1));
            fprintf(f,'%g\n',data(m,end));
        end

    fclose(f);
   


