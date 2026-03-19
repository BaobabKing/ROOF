%% Transcsipt composite_display


%% display composite image
function [Irgb,alpha] = composite_display(data,disp_gain)

Im = uint16(data);
L_H = zeros(2,4); % Initialize optimal image display range per image
TOL=[0.001 0.999]; %[0 1]; %; %[0 1]; %
for i=1:size(data,3)
    L_H(:,i) = stretchlim(Im(:,:,i),TOL);
    L_H(2,i) = max(L_H(2,i),0.0155);  % normalized signal intensity must be greater than 1000/(2^16-1)=0.0153
end

enhIm=Im;
for i=1:4
   tmp =  double(imadjust(Im(:,:,i),L_H(:,i),[]));
   if i==1 || i==2
      sortedtmp = sort(tmp(:),'descend');
      rng = sortedtmp(round(0.85*length(sortedtmp)));
   else
      rng = median(tmp(:)) + 3*mad(tmp(:));
   end
   enhIm(:,:,i)=(tmp>rng).*tmp +(tmp<=rng).*(rng) -rng;

end    
     
    
Irgb = cat(3, disp_gain(3)*enhIm(:,:,3),disp_gain(4)*enhIm(:,:,4),disp_gain(1)*enhIm(:,:,1)); % CK, Markers, Nucleus    
alpha = disp_gain(2)*enhIm(:,:,2);   % False
white =ones(size(alpha));
white =repmat(white,[1 1 3]);  % the two lines are equal to the below single line

figure; 
imshow(Irgb, 'InitialMag', 'fit'); axis on; impixelinfo;
truesize; hold on 
h = imshow(white); truesize; axis on;
set(h, 'AlphaData', alpha);   
hold off;           

%
% figure, 
% subplot(411); imshow(cat(3, enhIm(:,:,1), 0*enhIm(:,:,1),0*enhIm(:,:,1))); %imcontrast(gca); 
% subplot(412); imshow(enhIm(:,:,2),[]); imcontrast(gca);
% subplot(413); imshow(cat(3, 0*enhIm(:,:,3), 0*enhIm(:,:,3),enhIm(:,:,3))); %imcontrast(gca);
% subplot(414); imshow(cat(3, 0*enhIm(:,:,4), enhIm(:,:,4), 0*enhIm(:,:,4)));% imcontrast(gca);






