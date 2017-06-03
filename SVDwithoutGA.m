close all;
clc;
%loading the original image to be watermarked
orig_img = imread('lena_256.bmp');
figure, imshow(orig_img); title ('Original Image');
[M,N] = size(orig_img);
nblock=M/8;
%loading the image to be watermarked in the original image
wat_img = imread('njit_logo.jpg');
wat_img = imresize(wat_img, [nblock,nblock]); % make it 32 X 32
level = graythresh(wat_img);
wat_img = im2bw(wat_img,level); % make Black and White
figure, imshow(wat_img); title('WaterMark Image');
step_size = input('Step Size = ');         % ask step size size
%   Ask the type of Attacks 
fprintf('   1 - Midian Filter\n');
fprintf('   2 - Resize the image (scaling four times and then returned to original with bicubicinterpolation)\n');
fprintf('   3 - Salt and paper noise\n');
fprintf('   4 - Low pass filter\n');
fprintf('   5 - Image Jpeg Compression with Quality Factor 40\n');
attack_type = input('Attack Type: ');
    Wimg = SVD_Watermarked(orig_img, wat_img, step_size);
%   Attack 
%   Apply medianfilter
    if(attack_type==1)
        attacked_Image = medfilt2(Wimg);
    end
%   resize the image. first scaling (four times) then go back to original with 
%   bi-cubic interpolation 
    if(attack_type==2)
        attacked_Image = imresize(Wimg,4,'nearest');
        attacked_Image = imresize(attacked_Image, [M N]);
    end
%salt and paper noise
    if(attack_type==3)
        attacked_Image = imnoise(Wimg,'salt & pepper',0.002);
    end
%   Apply low pass filter
    if(attack_type==4)
        attacked_Image = uint8(conv2(double(Wimg), double(ones(3,3))/9));
        attacked_Image = imresize(attacked_Image, [M N]);
    end
    if(attack_type == 5)
        imwrite(Wimg,'Watermarked','JPEG','Quality',40);
        attacked_Image = imread('Watermarked.jpg');
    end
 % calculation of image quality degradiation after attack and inserting watermark
    [m,n] = size(orig_img);
    error = orig_img - attacked_Image;
    MSE = sum(sum(error^2))/(m*n);
    if (MSE > 0)
       peaksnr = 10*log10(255^2/MSE);
    else
       peaksnr=99;
    end
%   fprintf('The Peak Signal to noise ratio: %f db\n',peaksnr);

% Now extract the watermark image from attacked watermarked image
    wimg = extract_SVD_Watermarked(attacked_Image,step_size);
%Now calculate the similirity (Normalized correlation NC)between
%original watermark image and watermark image extracted from watermarked
%image 
%find out normalzed corelation
       norm_cor = corr2(wat_img, wimg);
%       m1 = mean2(orig_wat_img);
%       m2 = mean2(wimg);
%       s1 = sum(sum((orig_wat_img - m1)*(wimg - m2)));
%       s2 = sqrt(sum(sum((orig_wat_img - m1)^2))* sum(sum((wimg - m2)^2)));
%       norm_cor= s1/s2;
%fprintf('Normalized Cross correlation between original and extracted: %f \n',norm_cor);  
fprintf('PEAKSNR = %f NC = %f \n',peaksnr,norm_cor);
figure, imshow(attacked_Image);title('Attacked Image with Watermark');
figure, imshow(wimg);title('Watermark from attacked Image');
