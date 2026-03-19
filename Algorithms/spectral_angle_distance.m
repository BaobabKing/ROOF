function aSAD = spectral_angle_distance(W,Wg)
% Compute average spectral angle distance between estimated spectrum W and
% ground truth spectrum Wg.

[m,n]=size(W);

SAD = size(1,m);


for k=1:m

    w_k = W(k,:);
    wg_k = Wg(k,:);
    w_k = w_k/norm(w_k);
    wg_k = wg_k/norm(wg_k);
    SAD(k) = acos(sum(w_k.*wg_k));

end

aSAD = mean(SAD);