function Tot_Crosstalk = Heat_map_crosstalk(data,Titl, save_HeatMap)
% Input:
%      data   -   image 
%      Titl   -   heat map title, a string of characters 
%  save_HeatMap  - choose to save HeatMap to MATLAB data if it is flaged as 'Yes'
[rows,cols,~] = size(data);
X = zeros(rows*cols, 4,'uint16');

for i=1:size(data,3)
    tmp = data(:,:,i);
    X(:,i) = tmp(:);
end

R_Heat_Map = corrcoef(double(X));
nan_ind = isnan(R_Heat_Map);
R_Heat_Map(nan_ind) = 0;
for r=1:size(data,3)
    R_Heat_Map(r,r)=1; 
end

if strcmp(save_HeatMap,'Yes')
   save([Titl,'.mat'], 'R_Heat_Map');
end
% Calculate total crosstalk: as the sum of the off-diagonal elements of
% Heatmap squared minus number of diagonal elements
Tot_Crosstalk = sum(abs(R_Heat_Map(:))) - size(data,3); 
%Tot_Crosstalk = sum(sum(sqrt(R_Heat_Map.^2))) - size(data,3); 

figure, imagesc(R_Heat_Map);  title(Titl);

colormap turbo; %(hsv(512)); 
%colorbar; 
h = colorbar;
h.Limits = [0 1];
xlabel('Lambda (um)');
ylabel('Channel')
%xticks([1 2 3 4]);  % DAPI - 460 um; WBC - 490 um;  CK - 550 um;  Customers - 650 um.
xticklabels({'', '460','', '490','', '550','', '650'});
yticklabels({'', 'DAPI','', 'WBC','', 'CK','', 'Customers'});
%xticks([460 490 550 650]);  % DAPI - 460 um; WBC - 490 um;  CK - 550 um;  Customers - 650 um. 
%yticks([1 2 3 4]);  % DAPI, WBC, CK , Customers        
