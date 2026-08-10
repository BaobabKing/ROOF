function [ RGB, mask ] = KmeanClustering(Im, nColors)

disp_seg = 0;
Debug = 0;

ab = double(Im); % 1024*1024*2
nrows = size(ab,1);
ncols = size(ab,2);
nchannels=size(ab,3);
ab = reshape(ab,nrows*ncols,nchannels); 

% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',3); 
%%Label Every Pixel in the Image Using the Results from KMEANS
pixel_labels = reshape(cluster_idx,nrows,ncols);
if Debug==1
  figure, imshow(pixel_labels,[]), axis on; impixelinfo; title('image labeled by cluster index');
end
%%Create Images that Segment the H&E Image by Color.
segmented_images = cell(1,nColors);
for k = 1:nColors
    segmented_images{k} = pixel_labels==k;
end
if Debug==1
    figure;
    ax1 = subplot(2,3,1); imshow(segmented_images{1}), title('objects in cluster 1');
    ax2 = subplot(2,3,2); imshow(segmented_images{2}), title('objects in cluster 2');
    ax3 = subplot(2,3,3); imshow(segmented_images{3}), title('objects in cluster 3');
    ax4 = subplot(2,3,4); imshow(segmented_images{4}), title('objects in cluster 4');
    ax5 = subplot(2,3,5); imshow(segmented_images{5}), title('objects in cluster 5');
    linkaxes([ax1 ax2 ax3 ax4 ax5],'x');
end
% figure, imshow(segmented_images{3}), title('objects in cluster 3');
% figure, imshow(segmented_images{4}), title('objects in cluster 4');
% figure, imshow(segmented_images{5}), title('objects in cluster 5');
% 
if Debug==1
    L  = segmented_images{1} + 2*segmented_images{2} + 3*segmented_images{3} + 4*segmented_images{4} + 5*segmented_images{5};
    rgb = label2rgb(L,'jet',[.5 .5 .5]);
    figure, imshow(rgb); title('K-means clustering');impixelinfo;
end

%% Assign the cluster to the Background if it has the lowest
%%intensity across all channels  
[~,ind1] = min(cluster_center(:,1));
[~,ind2] = min(cluster_center(:,2));
[~,ind3] = min(cluster_center(:,3));
[~,ind4] = min(cluster_center(:,4));
G_ind = mode([ind1 ind2 ind3 ind4]);
background_seg = segmented_images{G_ind};
if disp_seg == 1
  figure, imshow(background_seg); axis on; title('Background');
end
% DAPI expression pixels corresponds to the row label of the highest 
% intensity in the first column of cluster_center.   
[~, N_ind ] = max(cluster_center(:,1));
Nucleus_seg = segmented_images{N_ind};
if disp_seg == 1
  figure, imshow(Nucleus_seg); axis on; title('Nucleus'); 
end

%% CK expression pixels corresponds to the row label of the highest 
% intensity element in the 3rd column of cluster_center.  
[~, CK_ind ] = max(cluster_center(:,3));
CK_seg = segmented_images{CK_ind};
if disp_seg == 1
   figure, imshow(CK_seg); axis on; title('CK'); 
end

%% WBC membrane expression pixels orresponds to the row label of the highest 
% intensity element in the 2nd column of cluster_center.  
[~, F_ind ] = max(cluster_center(:,2));
WBC_seg = segmented_images{F_ind};
if disp_seg == 1
  figure, imshow(WBC_seg); axis on; title('WBC membran marker'); 
end

%% Custom marker expression pixels orresponds to the row label of the highest 
% intensity element in the 4th column of cluster_center.  
[~,C_ind ] = max(cluster_center(:,4));
Custom_seg = segmented_images{C_ind};
if disp_seg == 1
   figure, imshow(Custom_seg); axis on; title('Customers'); 
end

%%% Pseudocolor the segments
%RGB = cat(3, 128*CK_seg, 128*Custom_seg, 128*Nucleus_seg );

RGB = cat(3, 128*CK_seg + 128*WBC_seg, 128*Custom_seg + 128*WBC_seg, 128*Nucleus_seg + 128*WBC_seg );
if Debug==1
  figure, imshow(RGB)
end

mask = cat(3, CK_seg, WBC_seg, Nucleus_seg, Custom_seg, background_seg);

end