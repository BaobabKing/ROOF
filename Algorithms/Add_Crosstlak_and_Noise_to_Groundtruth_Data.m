% This script was used to generate data contaiminated with crosstalk and noise 

cd ..;   s = cd; s0=s;  s = [s, '\Algorithms'];  cd Algorithms\

pathname_inpt = [s0,'\TestImages\Simulated\'];

whos('-file', [pathname_inpt,'simulated_ch_noise.mat']);

% Setup groundthruth channel crosstalk matrix (arranged in channel order of channel 1, channel 2, channel 3 and channel 4) (or DAPI, WBC, CK, Customer).
Wg = [1  0.005  0   0;   0.205  1  0.006 0.006; 0.0100 0.2750 1  0.0720; 0.1900 0.0120  0.0230  1];


load([pathname_inpt,'simulated_ch_noise.mat'], 'Ch_noises');  % To resuce randombessin results, use saved noise data 


%%% Load ground truth data
FileName = 'data_groundtruth.bin';
height=512; 
width=512; 
num_Channels=4;
imgSize=2*height*width;
header=0;
tail=0;
data_groundtruth = uint16(zeros(height,width,num_Channels));

fb = fopen([pathname_input filesep FileName],'rb');
for n=1:num_Channels
     % read in images
     fseek(fb, header + (n - 1)*(imgSize + tail),'bof');
     tmp = fread(fb, [width,height], 'uint16');
     data_groundtruth(:,:,n) = tmp;
end
fclose(fb);

%%% To test additional Poisson noise, uncomment the following 3 lines (jj for loop)
% for jj=1:num_Channels
%     data_groundtruth(:,:,jj) = imnoise(data_groundtruth(:,:,jj),"poisson");
% end

offset = 1000;

% Mix groundtruth data and add noise to simulate raw imaged data
data_raw = data_groundtruth; % initiallize raw data 
for k=1: num_channels

   tmp = Wg(k,1)*double(data_groundtruth(:,:,1)) + Wg(k,2)*double(data_groundtruth(:,:,2)) +...
         Wg(k,3)*double(data_groundtruth(:,:,3)) + Wg(k,4)*double(data_groundtruth(:,:,4));
               

   tmp1 = tmp + offset + Ch_noises(:,:,k);
   data_raw(:,:,k) = tmp1;

end

Pathname_save = pathname_input;
Filename_save = 'data_raw_ch1';
imwrite(uint16(data_raw(:,:,1)),[Pathname_save,Filename_save],'png');

Filename_save = 'data_raw_ch2';
imwrite(uint16(data_raw(:,:,2)),[Pathname_save,Filename_save],'png');

Filename_save = 'data_raw_ch3';
imwrite(uint16(data_raw(:,:,3)),[Pathname_save,Filename_save],'png');

Filename_save = 'data_raw_ch4';
imwrite(uint16(data_raw(:,:,4)),[Pathname_save,Filename_save],'png');


% Open a binary file for writing
Pathname_save = pathname_input;
fileID = fopen([Pathname_save filesep 'data_raw.bin'], 'wb');
% Write the matrix to the file
fwrite(fileID, data_raw, 'uint16'); % Use 'uint16' data type
fclose(fileID);
