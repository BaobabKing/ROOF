% EXAMPLE2     test run on real image.
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
CC_Methods ={'ROO', 'PICASSO', 'BINGO'};
Crosstalk_correction_method = CC_Methods{1}; %1 - 'ROOF'; % 2- 'PICASSO'; % 3- 'BINGO'.        

Segmented_area_only ='No'; 
if strcmp(Crosstalk_correction_method,'ROOF')
   Segmented_area_only = 'Yes';  
end


pathname_input = [s0,'\TestImages\Real\'];
Pathname_save = [ s0, '\ExperimentResults\Real\'];
save_corrected_data = 'Yes';
save_HeatMap = 'Yes';  
Save_data = 0; % '1' for 'To save', '0' for 'Not to save'

% Test image info  
num_Channels = 4;
header=0; tail=0;
height=401; width=351;
imgSize=2*height*width;
data_raw = uint16(zeros(height,width,num_Channels));

FileName = '40x_raw.dat'; 

%%% Load real data 
fb = fopen([pathname_input filesep FileName],'rb');
for n=1:num_Channels
     % read in images
     fseek(fb, header + (n - 1)*(imgSize + tail),'bof');
     tmp = fread(fb, [width,height], 'uint16');
     data_raw(:,:,n) = tmp';
end
fclose(fb);

decorrstretched_view_Dx(data_raw); title('data composite'); axis on;

% View channel images
figure('Name','raw CK'), imshow(data_raw(:,:,1),[]); title('DAPI'); imcontrast(gca); impixelinfo;
figure('Name','raw WBC'), imshow(data_raw(:,:,2),[]); title('WBC marker'); imcontrast(gca); impixelinfo;
figure('Name','raw DAPI'), imshow(data_raw(:,:,3),[]); title('CK'); imcontrast(gca); impixelinfo;
figure('Name','raw Customers'), imshow(data_raw(:,:,4),[]); title('Customers'); imcontrast(gca); impixelinfo;


%%% Raw data preprocessing: 
 % Estimate background and foreground mask 

[IRGB, mask]  = KmeanClustering(data_raw, 5);

Background_mask = mask(:,:,5);
Background = zeros(4,1);
for k = 1:num_Channels
    tmp = double(squeeze(data_raw(:,:,k)));
    Background(k) =  mean2(tmp(Background_mask)) + mad(tmp(Background_mask));

end
%%%% Alternatively, background mask can be segmented by thresholding
        % data_sum = sum(data_raw_cropped,3); 
        % th  = multithresh(data_sum,4);
        % Background_mask = data_sum<th(1);
        % Background_mask = imdilate(Background_mask,strel('disk',5));

% Subtract background from raw data
data_raw_bgrem = data_raw; % Initialize badk removed data
for i = 1:num_Channels
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

    X = zeros(height*width, num_Channels,'uint16');

    for i=1:num_Channels
        tmp = data_raw_bgrem(:,:,i);
        X(:,i) = double(tmp(:));
    end

end



% Setup crosstalk matrix from conventional routine calibration     
S = [1.0010  0.2069  0.0809  0.1945; 
     0.0050  1.0028  0.2772  0.0194; 
     0.0000  0.0062  1.0034  0.0232;
     0.0000  0.0065  0.0739  1.0018];

switch Crosstalk_correction_method 

    case 'ROOF'

        D = diag(S);
        E = [D D D D];
        U = S./E';
        Y = U';
        
        % Setup ROOF parameters
        alpha = 0.25;
        eps = 0.01; 
        Imax=2^16-1;
        maxiter = 15; 
        tStart = tic;
        [Xc,Wgc] = ROOF_unmixing(X, Y, alpha, eps, maxiter, size(data_raw_bgrem,1), size(data_raw_bgrem,2)); % adpative alpha 
         tElapsed = toc(tStart);
        
        
        % Convert flatten data back to images  
        data_raw_unmixed = zeros(height,width,num_Channels,'uint16');
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

            data_raw_unmixed = reshape(uint16(Xc),height, width, num_Channels);

    otherwise

        disp('Unkown method.');
        return;

end
% Visulize the corrected image by channels
figure('name','Unmixed Nucleus'), imshow(data_raw_unmixed(:,:,1),[]); axis on; impixelinfo; imcontrast(gca);
figure('name','Unmixed WBC marker'), imshow(data_raw_unmixed(:,:,2),[]); axis on; impixelinfo; imcontrast(gca); 
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
    

decorrstretched_view_Dx(data_raw_unmixed); title('Unmixed composite decor');


%%% Performonance quantification
Tot_Crosstalk_corrected = Heat_map_crosstalk(data_raw_unmixed,Titl, save_HeatMap);

Crosstalk_Correction_ratio = (Tot_Crosstalk_raw-Tot_Crosstalk_corrected)/Tot_Crosstalk_raw;
disp(['CC_ratio = ', num2str(Crosstalk_Correction_ratio)]);

% Compare the crosstalk Wgc from ROOF and the calibration crosstalk Y using
% control. NOte: Picasso and BINGO did not have calibration crosstalk to
% compare with. 
% 
if strcmp(Crosstalk_correction_method,'ROOF')
   aSAD = spectral_angle_distance(Wgc,Y);
   disp(['aSAD = ', num2str(aSAD)]);
end

disp(['Runtime (sec) = ', num2str(tElapsed)]) 

