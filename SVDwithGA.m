close all;
clc;
%loading the original image to be watermarked
orig_img = imread('lena_256.bmp');
figure, imshow(orig_img); title ('Original Image')
[M,N] = size(orig_img);
nblock=M/8;
%loading the image to be watermarked in the original image
wat_img = imread('njit_logo.jpg');
wat_img = imresize(wat_img, [nblock,nblock]); % make it 32 X 32
level = graythresh(wat_img);
wat_img = im2bw(wat_img,level); % make Black and White
figure, imshow(wat_img);  title('WaterMark Image');
% Now call genetic algo with 20 poulation, minimum value is 10, maximum
% value 90 with linear fitness function
npop = input('No of Iteration (Population Size) = ');         % ask population size
% with minimum value 10 and maximum value 90
[a,b,c] = genetic_algorithm (@(x)x^2,npop,npop,npop,npop,npop,npop,npop,1,150,250,.01);

%   Ask the type of Attacks 
fprintf('   1 - Midian Filter\n');
fprintf('   2 - Resize the image (scaling four times and then returned to original with bicubicinterpolation)\n');
fprintf('   2 - Salt and paper noise\n');
fprintf('   4 - Low pass filter\n');
fprintf('   5 - Image Jpeg Compression with Quality Factor 40\n');
attack_type = input('Attack Type: ');
max = 0.0; % max value of objective function
peaksnr_value = 0; % Initialize PEAKSNR value
NC = 0;  % Initialize Normalized Corelation value
step=0;
final_image = zeros(M,N);
for i=2:npop
    step_size = a(i); 
    Wimg = SVD_Watermarked(orig_img, wat_img, step_size);
% calculation of image quality degradiation after inserting watermark
    [m,n] = size(orig_img);
    error = orig_img - Wimg;
    MSE = (sum(sum(error^2)))/(m*n);
    if (MSE > 0)
       peaksnr = 10*log10(255^2/MSE);
    else
       peaksnr=99;
    end
%   fprintf('The Peak Signal to noise ratio: %f db\n',peaksnr);

%   Attack watermarked image
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
    if(attack_type==4)
        attacked_Image = uint8(conv2(double(Wimg), double(ones(3,3))/9));
        attacked_Image = imresize(attacked_Image, [M N]);
    end
    if(attack_type == 5)
        imwrite(Wimg,'Watermarked','JPEG','Quality',40);
        attacked_Image = imread('Watermarked.jpg');
    end
% Now extract the watermark image from attacked watermarked image
    wimg = extract_SVD_Watermarked(attacked_Image,step_size);
%Now calculate the similirity (Normalized correlation NC)between
%original watermark image and watermark image extracted from watermarked
%image 
       orig_wat_img = wat_img;
%find out normalzed corelation
       norm_cor = corr2(orig_wat_img, wimg);
%       m1 = mean2(orig_wat_img);
%       m2 = mean2(wimg);
%       s1 = sum(sum((orig_wat_img - m1)*(wimg - m2)));
%       s2 = sqrt(sum(sum((orig_wat_img - m1)^2))* sum(sum((wimg - m2)^2)));
%       norm_cor= s1/s2;
%fprintf('Normalized Cross correlation between original and extracted: %f \n',norm_cor);
% Calculate the objective function
       obfunc = norm_cor;
       fprintf('obfunc = %f step_size = %f PSNR= %f NC= %f\n',obfunc,step_size,peaksnr,norm_cor);
       if(obfunc > max)
           max=obfunc;
           step=step_size;
           peaksnr_value=peaksnr;
           NC = norm_cor;
           final_image=attacked_Image; % store the image into final_image
       end
end    
fprintf('Peaksnr = %f NC = %f \n',peaksnr_value,NC);
figure, imshow(final_image); title('Attacked Image with Watermark');
wimg = extract_SVD_Watermarked(final_image,step);
figure, imshow(wimg);title('Watermark from attacked Image');
