% This script was used to generate simulated data

cd ..;   s = cd; s0=s; 
pathname_input = [s0,'\TestImages\Simulated\'];


        Imax = 2^16-1;
        Pathname_save = pathname_input;
        Save_data = 1;
        
        h = 512;
        w = 512;
        num_channels = 4;
        data_groundtruth = zeros(h,w,num_channels, 'uint16');

        % Channel object MFI values
        MFI = [1000 400 200 600];
        SNR = [30 30 30 30];  % signal-to-noise ratio in each channel
        offset = 1000;

        % synthesize a CTC's nucleus channel image (groundtruth)        
        t=0:0.001:2*pi;
        x0=w/2;y0=h/2;
        R1x = 60; R1y=60;
        x = x0 + R1x*cos(t);
        y = y0 + R1y*sin(t);
        BW1 = poly2mask(x,y, h, w);
        data_groundtruth(:,:,1) = MFI(1)*double(BW1);   % ch1: DAPI
        
        % synthesize a CTC's CK channel image (ring shaped) groundtruth        
        R3x =110; R3y =80;
        x = x0+10 + R3x*cos(t);
        y = y0 + R3y*sin(t);
        BW3 = poly2mask(x,y, h, w);
        data_groundtruth(:,:,3) = MFI(3)*double(BW3-BW1); % Ch3: CK
%        BW1sm = imerode(BW1, strel('disk',1));
%        data_groundtruth(:,:,1) = MFI(1)*double(BW3-BW1sm);

        % synthesize a CTC's WBC channel image (ring shaped) groundtruth  
        R2x =120; R2y =90;
        x = x0 + R2x*cos(t);
        y = y0 + R2y*sin(t);
        BW2 = poly2mask(x,y, h, w);
        data_groundtruth(:,:,2) = MFI(2)*double(BW2-BW3);
        % BW3sm = imerode(BW3, strel('disk',1));
        % data_groundtruth(:,:,2) = MFI(2)*double(BW2-BW3sm);


        % synthesize a CTC's Customers channel image (ring shaped) groundtruth  
        BW4 = false(h,w);
%        rws = [220 250]; cls = [250 280];
        rws = [120 150]; cls = [250 280];
        BW4(rws,cls) = 1;
        BW4 = imdilate(BW4,strel('disk',10));
        data_groundtruth(:,:,4) = MFI(4)*double(BW4);

        Isum = data_groundtruth(:,:,1)+data_groundtruth(:,:,2)+data_groundtruth(:,:,3)+data_groundtruth(:,:,4);
        figure, imshow(Isum,[]); axis on; imcontrast(gca); impixelinfo;
        

       
        % Filename_save = 'groundtruth_data_ch1';
        % imwrite(uint16(data_groundtruth(:,:,1)),[Pathname_save,Filename_save],'png');
        % Filename_save = 'groundtruth_data_ch2';
        % imwrite(uint16(data_groundtruth(:,:,2)),[Pathname_save,Filename_save],'png');
        % Filename_save = 'groundtruth_data_ch3';
        % imwrite(uint16(data_groundtruth(:,:,3)),[Pathname_save,Filename_save],'png');
        % Filename_save = 'groundtruth_data_ch4';
        % imwrite(uint16(data_groundtruth(:,:,4)),[Pathname_save,Filename_save],'png');
      

       
        % Open a binary file for writing
        fileID = fopen([Pathname_save filesep 'data_groundtruth.bin'], 'wb');
        % Write the matrix to the file
        fwrite(fileID, data_groundtruth, 'uint16'); % Use 'uint16' data type
        fclose(fileID);
       