% EXAMPLE     test run on simulated image.
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

cd ..;   s = cd; s0=s;  s = [s, '\Algorithms']; path(s, path); cd Examples;
  
   
% Please select an algorithm ID
CC_Methods ={'ROOF','PICASSO', 'BINGO'};
Crosstalk_correction_method = CC_Methods{1}; % 1 - 'ROOF'; % 2- 'PICASSO'; % 3- 'BINGO'.        

Segmented_area_only ='No'; %'Yes'; % 
if strcmp(Crosstalk_correction_method,'ROOF')
   Segmented_area_only = 'Yes'; % 
end


pathname_input = [s0,'\TestImages\Simulated\'];
Pathname_save = [ s0, '\ExperimentResults\Simulated\'];
save_corrected_data = 'No'; 'Yes';
save_HeatMap = 'No';  

% image size 
Ht = 512; Wd = 512; % image size
num_channels = 4;


%%% Load smimulated data with crosstalk and noise
data_raw = zeros(Ht,Wd, num_channels,'uint16');
data_raw(:,:,1) = imread([pathname_input,'data_raw_ch1']);
data_raw(:,:,2) = imread([pathname_input,'data_raw_ch2']);
data_raw(:,:,3) = imread([pathname_input,'data_raw_ch3']);
data_raw(:,:,4) = imread([pathname_input,'data_raw_ch4']);


%%% Raw data preprocessing: 
 % Estimate background and foreground mask 

[IRGB, mask]  = KmeanClustering(data_raw, 5);

Background_mask = mask(:,:,5);
Background = zeros(4,1);
for k = 1:num_channels
    tmp = double(squeeze(data_raw(:,:,k)));
    Background(k) =  mean2(tmp(Background_mask)) + mad(tmp(Background_mask));

end
% Subtract background from raw data
data_raw_bgrem = data_raw; % Initialize badk removed data
for i = 1:num_channels
    data_raw_bgrem(:,:,i) = double(data_raw(:,:,i)) - Background(i);
end

        
Titl = 'Crosstalk heat map of raw image';
Tot_Crosstalk_raw = Heat_map_crosstalk(data_raw_bgrem,Titl, save_HeatMap);


if strcmp(Crosstalk_correction_method,'ROOF')

    %%% ROOF algorithm starts 
    % Flatten each band and create data matrix
    [rows, cols] = find(~Background_mask==1);
    num_pixels_foreground = numel(rows);
    I1 = zeros(num_pixels_foreground,1,'uint16');
    I2 = zeros(num_pixels_foreground,1,'uint16');
    I3 = zeros(num_pixels_foreground,1,'uint16');
    I4 = zeros(num_pixels_foreground,1,'uint16');
    
    for q=1:num_pixels_foreground
       I1(q) = data_raw_bgrem(rows(q),cols(q),1);
       I2(q) = data_raw_bgrem(rows(q),cols(q),2);
       I3(q) = data_raw_bgrem(rows(q),cols(q),3);
       I4(q) = data_raw_bgrem(rows(q),cols(q),4);
    end
    
    X = [I1 I2 I3 I4];

elseif strcmp(Crosstalk_correction_method,'BINGO')

    X = zeros(Ht*Wd, num_channels,'uint16');

    for i=1:num_channels
        tmp = data_raw_bgrem(:,:,i);
        X(:,i) = double(tmp(:));
    end

end

% Setup groundtruth channel crosstalk matrix (arranged in channel order: ch1,ch2,ch3 ,ch4. 
Wg = [1  0.005  0   0;   0.205  1  0.006 0.006; 0.0100 0.2750 1  0.0720; 0.1900 0.0120  0.0230  1];

switch Crosstalk_correction_method 

    case 'ROOF'

        % Puprposely modify the control estimated crosstalk to be different to groundtruth Wg.
        % Y = 0.75*Wg; % Undercorrection case
        % Y = Wg; % Ideal correction case
        Y = 1.2*Wg; % Overcorrection case
        for p=1:num_channels
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
        
        
        % Convert flatten data back to images  
        data_raw_unmixed = zeros(Ht,Wd,num_channels,'uint16');
        for q=1:num_pixels_foreground
            data_raw_unmixed(rows(q),cols(q),1) = Xc(q,1);
            data_raw_unmixed(rows(q),cols(q),2) = Xc(q,2);
            data_raw_unmixed(rows(q),cols(q),3) = Xc(q,3);
            data_raw_unmixed(rows(q),cols(q),4) = Xc(q,4);
                        
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

            data_raw_unmixed = reshape(uint16(Xc),Ht, Wd, 4);

    otherwise

        disp('Unkown method.');
        return;

end
% Visulize the corrected image by channels
figure('name','Unmixed Nucleus'), imshow(data_raw_unmixed(:,:,1),[]); axis on; impixelinfo; imcontrast(gca);
figure('name','Unmixed False'), imshow(data_raw_unmixed(:,:,2),[]); axis on; impixelinfo; imcontrast(gca); 
figure('name','Unmixed CK'), imshow(data_raw_unmixed(:,:,3),[]); axis on; impixelinfo; imcontrast(gca);
figure('name','Unmixed customers'), imshow(data_raw_unmixed(:,:,4),[]); axis on; impixelinfo; imcontrast(gca);

if strcmp(save_corrected_data,'Yes')

    Filename_save = [Crosstalk_correction_method,'_data_ch1'];
    imwrite(uint16(data_raw_unmixed(:,:,1)),[Pathname_save,Filename_save],'png');
    Filename_save = [Crosstalk_correction_method,'_data_ch2'];
    imwrite(uint16(data_raw_unmixed(:,:,2)),[Pathname_save,Filename_save],'png');
    Filename_save = [Crosstalk_correction_method,'_data_ch3'];
    imwrite(uint16(data_raw_unmixed(:,:,3)),[Pathname_save,Filename_save],'png');
    Filename_save = [Crosstalk_correction_method,'_data_ch4'];
    imwrite(uint16(data_raw_unmixed(:,:,4)),[Pathname_save,Filename_save],'png');

end
    
composite_display(uint16(data_raw_unmixed), [1 1 1 1]); title('Unmixed data composite'); axis on;

decorrstretched_view_Dx(data_raw_unmixed); title('Unmixed composite decor');


%%% Performonance quantification
% To compare RMSE, load groundthruth data
data_groundtruth = zeros(Ht,Wd, num_channels,'uint16');
data_groundtruth(:,:,1) = imread([pathname_input,'groundtruth_data_ch1']);
data_groundtruth(:,:,2) = imread([pathname_input,'groundtruth_data_ch2']);
data_groundtruth(:,:,3) = imread([pathname_input,'groundtruth_data_ch3']);
data_groundtruth(:,:,4) = imread([pathname_input,'groundtruth_data_ch4']);

Tot_Crosstalk_corrected = Heat_map_crosstalk(data_raw_unmixed,Titl, save_HeatMap);


SSIM = ssim(data_raw_unmixed,data_groundtruth);  
disp(['SSIM = ', num2str(SSIM)]);

aSAD = spectral_angle_distance(Wgc,Wg);
disp(['aSAD = ', num2str(aSAD)]);

aRMSE = averagedRMSE(data_raw_unmixed,data_groundtruth);
disp(['aRMSE = ', num2str(aRMSE)]);


disp(['Runtime (sec) = ', num2str(tElapsed)]) 

