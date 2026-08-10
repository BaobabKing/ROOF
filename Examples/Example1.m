% EXAMPLE:    Test on simulated image.
%
% 
% 
%
%   Author: Yongjian Yu 
%   Copyright (c) 2026 
%**************************************************
  
   
clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
workspace;	% Make sure the workspace panel is showing.

cd ..;   s = cd; s0 = s;  s = [s, '\Algorithms']; path(s, path); cd Examples;
  
   
% Please select an algorithm ID
CC_Methods ={'ROOF','PICASSO', 'BINGO'};
Crosstalk_correction_method = CC_Methods{2}; % 1 - 'ROOF'; % 2- 'PICASSO'; % 3- 'BINGO'.        
Y_factor = 1.2; % select from {0.75 1 1.2}; % To be used only for ROOF simulation experiments.




pathname_input = [s0,'\TestImages\Simulated\'];
Pathname_save = [ s0, '\ExperimentResults\Simulated\'];
save_corrected_data = 'No'; % 'Yes';% 
save_corrected_data_method = 'bin'; %'png'; %   

% image size 
height = 512; width = 512; % image size
num_Channels = 4;


%%% Load smimulated data with crosstalk and noise at different levels
SNR = 5; % or any of {30 25 20 15 10 5};
FileName = ['raw_data_SNR',num2str(SNR),'.bin'];
disp(FileName)

imgSize=2*height*width;
header=0;
tail=0;
data_raw = uint16(zeros(height,width,num_Channels));

fb = fopen([pathname_input filesep FileName],'rb');
for n=1:num_Channels
     % read in images
     fseek(fb, header + (n - 1)*(imgSize + tail),'bof');
     tmp = fread(fb, [width,height], 'uint16');
     data_raw(:,:,n) = tmp;
end
fclose(fb);

%%% Load groundtruth data
data_groundtruth = uint16(zeros(height,width,num_Channels)); 
fileName = 'groundtruth_data.bin';
fid = fopen([pathname_input filesep fileName],'rb');
for n=1:num_Channels
     % read in images
     fseek(fid, header + (n - 1)*(imgSize + tail),'bof');
     tmp = fread(fid, [width,height], 'uint16');
     data_groundtruth(:,:,n) = tmp;
end
fclose(fid);

%  
%%% For ROOF and BINGO algorithms to run properly, the channel index of the "for-correction" data needs to be in wavelength ascending order.
% The simulated data does not follow the right order requiring swapping of channels 1 and 3.
%%% Swap channels 1 and 3 
tmp = data_raw(:,:,1);
data_raw(:,:,1) = data_raw(:,:,3);
data_raw(:,:,3) = tmp;

tmp = data_groundtruth(:,:,1);
data_groundtruth(:,:,1) = data_groundtruth(:,:,3);
data_groundtruth(:,:,3) = tmp;

%%% Side note: PICASSO algorithm is not invariant with respect to channel ordering. 
% You may test to see: how Metric aSAR varies if channels 1 and 3 are swapped.


%%% Raw data preprocessing: 
 % Find background mask  
% Convert data_raw to MIP (Maximum Intensity Projection) 
I = max(data_groundtruth,[],3);
Background_mask = I==0;
% figure, imshow(Background_mask);


% Estimate background for background removal
Background = zeros(4,1);
for k = 1:num_Channels
    tmp = double(squeeze(data_raw(:,:,k)));
    Background(k) =  mean2(tmp(Background_mask)) + mad(tmp(Background_mask));

end
% Subtract background from raw data (resulting in "for-correction data").
data_raw_bgrem = data_raw; % Initialize badk removed data
for i = 1:num_Channels
    data_raw_bgrem(:,:,i) = double(data_raw(:,:,i)) - Background(i);
end

 

% Visulize the for-correction image by channels
figure('name','Nucleus'), imshow(data_raw_bgrem(:,:,1),[]); axis on; impixelinfo; imcontrast(gca);
figure('name','WBC marker'), imshow(data_raw_bgrem(:,:,2),[]); axis on; impixelinfo; imcontrast(gca); 
figure('name','CK'), imshow(data_raw_bgrem(:,:,3),[]); axis on; impixelinfo; imcontrast(gca);
figure('name','Customers'), imshow(data_raw_bgrem(:,:,4),[]); axis on; impixelinfo; imcontrast(gca);


if strcmp(Crosstalk_correction_method,'ROOF')

    % Reshape 3D data into 2D by flattening each channel into a column vector
    X0 = reshape(data_raw_bgrem, [], num_Channels);
        
    % Extract only foreground pixels
    [rows, cols] = find(Background_mask==0);
    num_pixels_foreground = numel(rows);
    Foreground_mask = ~Background_mask(:);
    X = zeros(num_pixels_foreground,num_Channels);

    for j=1:num_Channels
       X(:,j) = X0(Foreground_mask,j);
    end

elseif strcmp(Crosstalk_correction_method,'BINGO')
    
    % Reshape 3D data into 2D by flattening each channel into a column vector
    X = reshape(data_raw_bgrem, [], num_Channels);

end

% Setup groundtruth channel crosstalk matrix (arranged in channel order: ch1, ch2,ch3 ,ch4. 
Wg = [1  0.005  0   0;   0.205  1  0.006 0.006; 0.0100 0.2750 1  0.0720; 0.1900 0.0120  0.0230  1];
% Wg =[1  0.006  0  0.023; 0.275 1 0.005 0.012;0.01 0.205 1 0.190; 0.072 0.006 0 1];
% Wg = crosstalk_refomulate_if_twochannels_swarpped(Wg', 1, 3);

switch Crosstalk_correction_method 

    case 'ROOF'

        % Puprposely modify the control estimated crosstalk to be different to groundtruth Wg.
        Y = Y_factor*Wg; % Y_factor = 0.75 Undercorrection case
                         % Y_factor = 1; Ideal correction case
                         % Y_factor = 1.2; Overcorrection case
        for p=1:num_Channels
           Y(p,p) = 1;  % set each diagonal element in Y to unity. 
        end
        
        % Setup ROOF parameters
        alpha = 0.25;
        eps = 0.01; 
        Imax=2^16-1;
        maxiter = 15; 
        tStart = tic;
        [Xc,Wgc] = ROOF_unmixing(X, Y, alpha, eps, maxiter, size(data_raw_bgrem,1), size(data_raw_bgrem,2)); % adpative alpha 
        tElapsed = toc(tStart);
        
        
        % Convert the corrected flattened data back to 3D image  
        data_raw_unmixed = zeros(height,width,num_Channels,'uint16');
        
        for kk=1:num_Channels
            for q=1:num_pixels_foreground
                data_raw_unmixed(rows(q),cols(q),kk) = Xc(q,kk);              
            end
        end

    case 'PICASSO'  % Junnyoung Seo, Yeonbo Sim et al 2022

            per_bg = 0; %% percentile value (0~100) of background
            step_size = 0.2; %% step_size of updating alpha
            maxIter = 25; %% the number of iteration

            
            qN = 100;
            tStart = tic;
            [data_raw_unmixed, unmixing_log] = PICASSO_4C(data_raw_bgrem, qN, maxIter, step_size, 0);
            Wgc = inv(unmixing_log.UMM);
            tElapsed = toc(tStart);
            % Reshape 2D corrected data back to 3D tensor form
            data_raw_unmixed = uint16(data_raw_unmixed);

    case 'BINGO'
            
            alpha = 0.25;
            eps = 0.01;
            maxiter = 5;
            mu = 0.01;   % gradient descent step size
            r = size(X,2);  % number of fluorophores
            spH = 0.5;  % Target sparseness, value range [0 1]
            tStart = tic;
            [Xc, Wgc] = BINGO_unmixing(X, r, alpha, mu, eps, maxiter, spH);
            tElapsed = toc(tStart);

            data_raw_unmixed = reshape(uint16(Xc),height, width, 4);

    otherwise

        disp('Unkown method.');
        return;

end

% Visulize the corrected image by channels
figure('name','Corrected Nucleus'), imshow(data_raw_unmixed(:,:,1),[]); axis on; impixelinfo; imcontrast(gca);
figure('name','Corrected WBC marker'), imshow(data_raw_unmixed(:,:,2),[]); axis on; impixelinfo; imcontrast(gca); 
figure('name','Corrected CK'), imshow(data_raw_unmixed(:,:,3),[]); axis on; impixelinfo; imcontrast(gca);
figure('name','Corrected Customer marker'), imshow(data_raw_unmixed(:,:,4),[]); axis on; impixelinfo; imcontrast(gca);

if strcmp(save_corrected_data,'Yes')
    %%% save as png or bin format
    if strcmp(save_corrected_data_method,'png')
        Filename_save = [Crosstalk_correction_method,'_data_ch1'];
        imwrite(uint16(data_raw_unmixed(:,:,1)),[Pathname_save,Filename_save],'png');
        Filename_save = [Crosstalk_correction_method,'_data_ch2'];
        imwrite(uint16(data_raw_unmixed(:,:,2)),[Pathname_save,Filename_save],'png');
        Filename_save = [Crosstalk_correction_method,'_data_ch3'];
        imwrite(uint16(data_raw_unmixed(:,:,3)),[Pathname_save,Filename_save],'png');
        Filename_save = [Crosstalk_correction_method,'_data_ch4'];
        imwrite(uint16(data_raw_unmixed(:,:,4)),[Pathname_save,Filename_save],'png');
    
    elseif strcmp(save_corrected_data_method,'bin')
        
        if strcmp(Crosstalk_correction_method,'ROOF')
            Filename_save = [Crosstalk_correction_method,'_corrected_data_SNR',num2str(SNR),'_Y',num2str(Y_factor), '.bin'];
        else
            Filename_save = [Crosstalk_correction_method,'_corrected_data_SNR',num2str(SNR(1)),'.bin'];
        end
        fid = fopen([Pathname_save, Filename_save],'w');
        fwrite(fid,data_raw_unmixed,'uint16');
        fclose(fid);
    end
end
    

%%% Performonance quantification
SSIM = ssim(data_raw_unmixed,data_groundtruth);  
disp(['SSIM = ', num2str(SSIM)]);

aSAD = spectral_angle_distance(Wgc,Wg);
disp(['aSAD = ', num2str(aSAD)]);

aRMSE = averagedRMSE(data_raw_unmixed,data_groundtruth);
disp(['aRMSE = ', num2str(aRMSE)]);


disp(['Runtime (sec) = ', num2str(tElapsed)]) 

% End of test