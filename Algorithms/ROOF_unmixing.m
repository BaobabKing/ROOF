function [Dc, W] = ROOF_unmixing(D, Y, alpha, eps, maxiter, rows, cols)
% The generic NMF Unmixing Algorithm 1.
%  Input:  
%        D   -   preprocessed image data matrix  dimension: rows*cols x num_channels
%        Y   -   reference crosstalk matrix
%    alpha   -   regularization weight
%      eps   -   convergence tolerance
%  maxiter   -   num. of maximaum iterations
%     rows   -   number of rows of raw channel image
%     cols   -   number of cols of raw channel image
%   Output   
%       Dc   -   crostalk corrected image data matrix
%        W   -   optimized corosstalk matrix
%   
% This version 2 makes the value of alpha adaptive. 03/10/2025 
%**********************************************************
debug = 0;

 if debug == 1
    D_unflattened = reshape(D,rows,cols,4);
    D_unflattened = uint16(D_unflattened);
 end

Imax = 2^16-1;

D = double(D)/Imax;
[m,n] = size(D);
k = n;

insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end));

% Initialize the W matrix .
W = Y; % W represents the spatial concentration distribution of each fluorophore.


% H - the feature matrix reflects the intensity distribution of each fluorophore in each 
% channel, corresponding to its estimated spectrum

sse_diff = 1;
sse = 10;
niter = 0;
%maxiter = 2;
small0 = eps^3;
small = 0;
alpha_r = zeros(1,n);
alpha_r(1:n) = alpha; 
% Iteration
while sse_diff > small0
    if niter > maxiter
        break
    end    
    niter = niter + 1;
    % Calculate F under nonnegativity constraint:
    if niter==1
        F = D/W;
        nan_ind = isnan(F);
        inf_ind = isinf(F);
        neg_ind = F<=0;
        F(nan_ind) = small;
        F(inf_ind) = small;
        F(neg_ind) = small;

    else
        F = F.*(D*W')./(F*W*W');
        nan_ind = isnan(F);
        inf_ind = isinf(F);
        neg_ind = F<=0;
        F(nan_ind) = small;
        F(inf_ind) = small;
        F(neg_ind) = small;
    end
    if debug == 1
        F_unflattened = reshape(F,rows,cols,4);
        F_unflattened = uint16(F_unflattened*Imax);
        % Swap channels 1 and 3 to alighn data coulmn order with the wavelength Index
        tmp = F_unflattened(:,:,3);
        F_unflattened(:,:,3) = F_unflattened(:,:,1);
        F_unflattened(:,:,1)= tmp;
        composite_display(F_unflattened, [1 1 1 1]); title('raw data composite'); axis on;
    end
    % nan_ind = isnan(F);
    % inf_ind = isinf(F);
    % neg_ind = F<=0;
    % F(nan_ind) = small;
    % F(inf_ind) = small;
    % F(neg_ind) = small;
    % Normalize the rows. 
    % for nrows = 1:k
    %     F(nrows,:) = F(nrows,:)/norm(F(nrows,:));
    % end
    % Here we check again for NaN, inf, or negative values before the
    % pseudoinverse calculation.
    % nan_ind = isnan(F);
    % inf_ind = isinf(F);
    % neg_ind = F<=0;
    % F(nan_ind) = small;
    % F(inf_ind) = small;
    % F(neg_ind) = small;

 
    % Update regularization paremeter alpha_r
    for r = 1:n
       F_minus_r = F;
       F_minus_r(:,r) =[];
       % find eigenvalues andeigenvectors of  (F_minus_r)' * F_minus_r
       [Q,V] = eig((F_minus_r)'* F_minus_r);
       beta = diag(V); 
       Dr = D(:,r);
       Fr = F(:,r);       
       Yr = Y(:,r);
       Yr(r)=[];

       Z = Q'*(F_minus_r)'*(Dr-Fr-F_minus_r*Yr);

       H1 = sum(Z.^2./(beta+alpha_r(r)));
       H2 = sum(Z.^2./(beta+alpha_r(r)).^2);

       alpha_r(r) = sum(beta./(beta + alpha_r(r)))*(H1/H2 + alpha_r(r));

    end
    alpha = mean(alpha_r);


    
    % Calculate W:
    W_tmp = [];
   
    for r = 1:n 

       Dr = D(:,r);
       Fr = F(:,r);       
       Yr = Y(:,r);
       Yr(r)=[];
       F_minus_r = F;
       F_minus_r(:,r) =[];
       
%       Wr1 = inv(F_minus_r'*F_minus_r + alpha_r(r)*eye(n-1))*(F_minus_r'*(Dr-Fr) + alpha_r(r)*eye(n-1)*Yr);
       Wr1 = inv(F_minus_r'*F_minus_r + alpha*eye(n-1))*(F_minus_r'*(Dr-Fr) + alpha*eye(n-1)*Yr);
       
       % Insert 1 back at rth row in Wc1
       Wr2 = insert(1, Wr1', r-1);
       W_tmp = [W_tmp Wr2'];
    end

    W = W_tmp;

    % Implement nonnegativity:
    nan_ind = isnan(W);
    inf_ind = isinf(W);
    neg_ind = W<=0;
        
    W(nan_ind) = small;        
    W(inf_ind) = small;        
    W(neg_ind) = small;   
    

    % Calculate residual:
    D_hat = F*W;
    R = D - D_hat;
    sse_old = sse;
    sse = R(:)'*R(:) + alpha*norm(W(:)-Y(:)); 
    disp(['niter=',num2str(niter)])
    sse_diff = abs(sse_old - sse)
    Dc = F*Imax;
end


