function [Wimg]= SVD_Watermarked(orig_img, wat_img, step_size)
%output 
% - wimg : watermarked image
%Input
% - orig_image : original image without watermarking
% - wat_img : watermark image
% - step_size : The interval value
[M,N] = size(orig_img);
orig_img = double(orig_img);
ind=(M*N)/64;
dlarge = zeros(ind,1); % create array of size ind;
block = zeros(8,8,ind); % create and initialize block
index=0;
for i = 1:8:M;
    for j = 1:8:N;
        index=index+1;
        block(:,:,index) = orig_img(i:i+7,j:j+7);  % make three dimensional array
        [Uimg,Simg,Vimg] = svd(block(:,:,index));
        dlarge(index)=Simg(1,1);              % store the the large value of S
    end
end
%max = max(dlarge);
max = dlarge(1);          % Find the maximum value from dlarge array
for jj = 2:numel(dlarge)
    if max < dlarge(jj)
        max = dlarge(jj);
    end
end
%min = min(dlarge);
min = dlarge(1);          % find the minimum value from dlarge array
for jj = 2:numel(dlarge)
    if min > dlarge(jj)        
        min = dlarge(jj);
    end
end
interval_table = min-step_size : step_size : max+step_size;      % construct the table 
[row col]=size(interval_table);
for i=1:index
    for j=1:col-1
        if ((dlarge(i)> interval_table(j)) && (dlarge(i) < interval_table(j+1)));
            if(wat_img(i)==1);   % change the dlarge value based on watermarking image
                dlarge(i)= (interval_table(j)+(interval_table(j)+interval_table(j+1))/2.)/2.;
            else
                dlarge(i)= (interval_table(j+1)+(interval_table(j)+interval_table(j+1))/2.)/2.;
            end
        end
    end
end
index=0;
Wimg = zeros(8,8,ind); % create and initialize Wimg
for i = 1:8:M;
    for j = 1:8:N;
        index=index+1;
        block(:,:,index) = orig_img(i:i+7,j:j+7);
        [Uimg,Simg,Vimg] = svd(block(:,:,index));
        Simg(1,1)= dlarge(index);       %put the large value based on watermarked image
        Wimg(:,:,index)=Uimg*Simg*Vimg'; %inverse the SVD
    end
end  
%Wimg is three dimensional matrix we have to change to two dimensional
%matrix
nrow = M/8;
ncol = N/8;
D = zeros(8,8,ind); % create and initialize D
final = zeros(8,8,ind); %create and initialize D
index=0;
for row=1:nrow
    index=index+1;
    D=Wimg(:,:,index);
    for col=2:ncol
       index=index+1;
       D=horzcat(D,Wimg(:,:,index));
%       fprintf('row= %d col= %d',row,col);
    end
    if(row==1)
        final=D;
    else
        final=vertcat(final,D);
    end
end
Wimg = uint8(final);


  

