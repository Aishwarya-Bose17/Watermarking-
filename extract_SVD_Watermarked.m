function [wat_img]= extract_SVD_Watermarked(watermarked_img,step_size)

%output 
% - wat_img : watermark image
%Input
% - watermarked_image : watermarked image
% - step_size : The interval value

[M,N] = size(watermarked_img);
watermarked_img = double(watermarked_img);
%dlarge = [];
%block=[];
nb=(M*N)/64;
dlarge = zeros(nb,1); % create array of size nb;
block = zeros(8,8,nb); % create and initialize block
ind=0;
for i = 1:8:M;
    for j = 1:8:N;
        ind=ind+1;
        block(:,:,ind) = watermarked_img(i:i+7,j:j+7);  % make three dimensional array
        [Uimg,Simg,Vimg] = svd(block(:,:,ind));
        dlarge(ind)=Simg(1,1);              % sdtore the the large value of S
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
interval_table = min-step_size : step_size : max+step_size;           % construct the table 
%disp(interval_table);
[row col]=size(interval_table);
wat_img = zeros(nb,1); % create array of size nb
for i=1:ind
    for j=1:col-1
        if ((dlarge(i)> interval_table(j)) && (dlarge(i) < (interval_table(j)+ interval_table(j+1))/2))
            wat_img(i)=1;   
        elseif ((dlarge(i)> (interval_table(j)+interval_table(j+1))/2) && (dlarge(i) < interval_table(j+1)))
            wat_img(i)=0;
        end
    end
end
nb=sqrt(ind);
% To convert 1D array to 2D array in major row
%for i=1:nb
%  for j=1:nb
%     final(i,j)=wat_img((j-1)*nb+i);
%  end
%end
wat_img = reshape(wat_img,nb,nb); % In matlab, function is there to do this. Above code does the same


  

