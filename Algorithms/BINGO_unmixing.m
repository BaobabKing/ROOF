function [Xc, H] = BINGO_unmixing(X, r, alpha, mu, eps, maxiter, spH)
% The sparse NMF Unmixing Algorithm 1.
% X~=W*H
% Input
%       X   - image matrix
%       r   - bumnber of fluorophores 
%   alpha   - sparse regularization weight
%      mu   - gradient descent step size
%     eps   - convergence check bound
% maxiter   - maximum number of iterations


Imax=2^16-1;
X = double(X)/Imax;
sse_diff = 1;
sse = 1;
niter = 0;

small0 = eps^3;
small = eps;

% Initalization using SVD 
[W0, H0] = BINGO_NNDSVD_Initialization(single(X),r);
nan_ind = isnan(W0);
inf_ind = isinf(W0);
neg_ind = W0<=0;
W0(nan_ind) = small;
W0(inf_ind) = small;
W0(neg_ind) = small;

% Apply projected gradient descent

% Iteration
while sse_diff > small0
    if niter > maxiter
        break
    end    
    niter = niter + 1;

    [W,H, sparseness, sparseness_r] = Project_gradient_descent_update(X, W0, H0, alpha, mu, spH);


    nan_ind = isnan(W);
    inf_ind = isinf(W);
    neg_ind = W<=0;
    W(nan_ind) = small;
    W(inf_ind) = small;
    W(neg_ind) = small;



    % Calculate residual:
    X_hat = W*H;
    R = X - X_hat;
    sse_old = sse;
    sse = R(:)'*R(:); 
    sse_diff = abs(sse_old - sse);
    W0 = W;
    H0 = H;

end

Xc = W*Imax;

% [m,n] = size(X);
% k = r;
% 
% % Initialize the W matrix and normalize the columns.
% W = rand(m,k); % W represents the spatial concentration distribution of each fluorophore.
% for ncols = 1:k
%     W(:,ncols) = W(:,ncols)/norm(W(:,ncols));
% end
% 
% % H - the feature matrix reflects the intensity distribution of each fluorophore in each 
% % channel, corresponding to its estimated spectrum
% 
% sse_diff = 1;
% sse = 1;
% niter = 0;
% maxiter = 10;
% small0 = eps;
% small = 0;
% % Iteration
% while sse_diff > small0
%     if niter > maxiter
%         break
%     end    
%     niter = niter + 1;
%     % Calculate H and implement nonnegativity:
%     H = pinv(W)*X;
%     nan_ind = isnan(H);
%     inf_ind = isinf(H);
%     neg_ind = H<=0;
%     H(nan_ind) = small;
%     H(inf_ind) = small;
%     H(neg_ind) = small;
%     % Normalize the rows. 
%     for nrows = 1:k
%         H(nrows,:) = H(nrows,:)/norm(H(nrows,:)); 
%         % H(nrows,nrows)=1;
%     end
%     % Here we check again for NaN, inf, or negative values before the
%     % pseudoinverse calculation.
%     nan_ind = isnan(H);
%     inf_ind = isinf(H);
%     neg_ind = H<=0;
%     H(nan_ind) = small;
%     H(inf_ind) = small;
%     H(neg_ind) = small;
% 
%     % Calculate A:
%     W = X*pinv(H);
% 
%     % Implement nonnegativity:
%     nan_ind = isnan(W);
%     inf_ind = isinf(W);
%     neg_ind = W<=0;
% 
%     W(nan_ind) = small;        
%     W(inf_ind) = small;        
%     W(neg_ind) = small;        
% 
% 
%     % Calculate residual:
%     X_hat = W*H;
%     R = X - X_hat;
%     sse_old = sse;
%     sse = R(:)'*R(:); 
%     sse_diff = abs(sse_old - sse);
% end


function [W0, H0] = BINGO_NNDSVD_Initialization(X,r)
% NNDSVD Initialization 
% Input：
%    X    - Input the image matrix X 
%    r    - the number of fluorophores to be solved.
% Output：
%    W0 and H0 (the initialization value W0 and H0).
%   
% Algorithm 2:
% At first, it performs a singular 
% value decomposition of the image matrix X with r as the rank, and the goal of the second SVD is to 
% calibrate all elements of the matrix to non-negative. 
%

% 1. Compute the largest r singular triplets of X: X = U*S*V'
% Singular value decomposition expresses an m-by-n matrix X as X = U*S*V'. 
% Here, S is an m-by-n diagonal matrix with singular values of X on its diagonal. 
% The columns of the m-by-m matrix U are the left singular vectors for corresponding singular values. 
% The columns of the n-by-n matrix V are the right singular vectors for corresponding singular values. 
% V' is the Hermitian transpose (the complex conjugate of the transpose) of V.

[m,n] = size(X);

[U,S,V] = svd(X, "econ");
W0 = zeros(m,r);
H0 = zeros(r, n);

for j=1:r

    u = U(:,j);
    v = V(:,j);
    u_plus = u.*(u>=0) + 0*u.*(u<0);   %  u_plus is written as all the elements of ≥ 0 in vector u, and the rest are replaced by 0.
    u_minus = u.*(u<0) + 0*u.*(u>=0);  % u_minus is written as all the elements of < 0 in vector u, and the rest are replaced by 0.
    v_plus = v.*(v>=0) + 0*v.*(v<0);
    v_minus = v.*(v<0) + 0*v.*(v>=0);
    
    if (norm(u_plus)*norm(v_plus)) > (norm(u_minus)*norm(v_minus))
        
        D = sqrt(S(j,j)*norm(u_plus)*norm(v_plus));
        W0(:,j) = D/norm(u_plus)*u_plus;
        H0(j,:) = D/norm(v_plus)*(v_plus)';


    else
        D = sqrt(S(j,j)*norm(u_minus)*norm(v_minus));
        W0(:,j) = D/norm(u_plus)*u_plus;
        H0(j,:) = D/norm(v_minus)*(v_minus)';

    end


end


function [W,H, sparseness, sparseness_r] = Project_gradient_descent_update(V, W0, H0, alpha, mu, spH)
% Input:
%      V    -  Data matrix   ~W*H
%     W0    -  concentration W at k-th iterate
%     H0    -  spectra H at k-th iterate 
%  alpha    -  Sparsity penality weighting factor. [0, 1]
%     mu    -  gradient descent step size
%    spH    -  spH is an adjustable parameter ranging from 0 to 1.
% Output:
%      W    -  One-iterate updated W
%      H    -  one-iterate updated H
%
% Algorithm 3 ：Project gradient descent update rule
% This combined step consists of a gradient descent step (Step 1) followed by 
% projection onto the closest point satisfying both the non-negativity and the 
% unit-norm constraints (Steps 2 and 3). Where .∗ and ./ denote elementwise 
% multiplication and division, respectively.
%***********************************************************


n = size(V,2);
V = double(V);
% Step 1. Set H1= H - mu*W'*(W*H-V) 
H1 = H0 - mu*W0'*(W0*H0-V);   % a n-by-n matrix
% Step 2. Project each row of H to be nonnegative, have unit L2 norm, L1 norm set 
% to achieve desired sparseness.
H1 = max(H1,0);
H = H1;
for k = 1:n  % loop through each fluorophore
       
       %H1(k,:) = H1(k,:)/norm(H1(k,:)); 
              
       L1 = sum(abs(H1(k,:)));    
       L2 = norm(H1(k,:)); % L2 equals 1
       if L2==0
          H1(k,k)=1;
       end
       sparseness =  (sqrt(n) - L1/L2)/(sqrt(n)-1); % Note that the spareser the vector H is, the larger such a sparseness is valued.  
        
       H2 = H1(k,:);
       % Divid H2 by a factor B (:= sqrt(n)-(sqrt(n)-1)*spH) and restore the maximum in H2 to max(H1(k,:))
       B = sqrt(n)-(sqrt(n)-1)*spH;
       B = B^alpha; % if alpha=0, B=1, there is no sparse regulizaiton.
       H2 = H2/B;
     
       % restore the maximum in H2 to max(H1(k,:))
       [~,ind] = max(H2);
       H2(ind) = max(H1(k,:));

       L1_r = sum(abs(H2));
       L2_r = norm(H2);
       sparseness_r =  (sqrt(n) - L1_r/L2_r)/(sqrt(n)-1);
       
       H(k,:) = H2;  
    
    

end

% 3. W := W*(V*H')./(W*H*H')   The multiplicative update rule by Daniel D. Lee
% and H.S. Seung, Algorithm for non-negative Matrix Factorization
W = W0.*(V*H')./(W0*(H*H'));