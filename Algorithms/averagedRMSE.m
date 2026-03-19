function aRMSE = averagedRMSE(I1,I2)
% Compute average spectral angle distance between estimated spectrum W and
% ground truth spectrum Wg.

[m,n,d]=size(I1);
S = m*n;

RMSEs = size(1,d);


for k=1:d

%   RMSEs(k) = sqrt(immse(I1(:,:,k),I2(:,:,k)));
   X1 = double(I1(:,:,k));
   X2 = double(I2(:,:,k));
   RMSEs(k) = sqrt(sum((X1(:)-X2(:)).^2)/S);

end

aRMSE = mean(RMSEs);