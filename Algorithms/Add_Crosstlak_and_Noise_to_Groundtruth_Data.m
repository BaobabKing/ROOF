% This script was used to generate data contaiminated with crosstalk and noise 


% Setup groundthruth channel crosstalk matrix (arranged in channel order of channel 1, channel 2, channel 3 and channel 4) (or DAPI, WBC, CK, Customer).
Wg = [1  0.005  0   0;   0.205  1  0.006 0.006; 0.0100 0.2750 1  0.0720; 0.1900 0.0120  0.0230  1];

load('simulated_ch_noise.mat', 'Ch_noise');  % To resuce randombessin results, use saved noise data 

% Mix groundtruth data and add noise to simulate raw imaged data
data_raw = data_groundtruth; % initiallize raw data 
for k=1: num_channels

   tmp = Wg(k,1)*double(data_groundtruth(:,:,1)) + Wg(k,2)*double(data_groundtruth(:,:,2)) +...
         Wg(k,3)*double(data_groundtruth(:,:,3)) + Wg(k,4)*double(data_groundtruth(:,:,4));
               

   tmp1 = tmp + offset + Ch_noises(:,:,k);
   data_raw(:,:,k) = tmp1;

end
