function [Kappa, Gamma, R_square, R_cv, W ,M] = mkmeans95(Map_gfp, K, iter, b, lamda, varargin)
%% Modified K-means clustering (Pascual-Marqui RD et.al. IEEE Trans Biomed Eng. 1995;42:658¨C65.)
%  Clustering algorithm and smoothing
% 
%  inputs:  1) Map_gfp :  Maps of GFP peaks (M*N matrix) 
%                         M - number of electordes  N - number of time frame
%           2) K       :  Number of Microstates
%           3) iter    :  itertions for optimum solution 
%           4) b       :  Section length
%           5) lamda   :  Penalty parameter for nonsmoothness
%           6) Default :  tol,iter_thereshold
%
%  outputs: 1) Kappa   : Series (Each peak belongs to its corresponding cluster)
%           2) Gamma   : Microstate (M*N matrix)
%                        M - number of electrodes N - number of micorstate
%           3) R_square: the extent of fitting (1.00 is the best)
%           4) R_cv    : Variance in cross validation
%           5) W       : KL criterion param (dispersion of clusters)
%           6) M       : KL criterion param
%
% Anthor: Allard Wen Shi  
% Copyright (C) Allard Shi from Xian Jiaotong University  

%% Default Params
if ~exist('tol')
    tol = 1e-6;
end

if ~exist('iter_thereshold')
    iter_thereshold = 10000;
end
%% Initialize the params
sigma_o = 0; sigma_u = 1000; sigma = 1000;
[num_channels,frame] = size(Map_gfp); 
Gamma = zeros(num_channels,K);
gamma = zeros(num_channels,K);
Nbkt = zeros(frame-2*b,K);
GEV = zeros(K,1);                              % Global Explained Variance for each microstate
GFP = std(Map_gfp,1);                          % Global Field Power
GFP_all = sum(GFP.^2);                         % Sum of GFP

for t = 1:iter
    for i = 1:K
        gamma(:,i) = randn(num_channels,1);                     
        gamma(:,i) = gamma(:,i)/norm(gamma(:,i),2);            % |gamma| =1
    end

    %% Circulation Optimization
    while abs(sigma_o-sigma_u) > sigma_u*tol
        sigma_o = sigma_u;
        [~,L] = max((Map_gfp'*gamma).^2,[],2);  
        for k = 1:K
            index = find(L==k);
            S = Map_gfp(:,index)*Map_gfp(:,index)';            % compute Sk
            [evector,evalue] = eig(S);
            [~,lam] = max(diag(evalue));
            gamma(:,k) = evector(:,lam);
            gamma(:,k) = gamma(:,k)/norm(gamma(:,k),2);        % Normalization
            GEV(k) = GFP(1,index).^2*compute_spatial_correlation(Map_gfp(:,index),gamma(:,k)).^2/GFP_all;
        end

        % Compute new sigma_u
        sigma_u = 0;
        for j = 1:frame
            sigma_u = sigma_u + (Map_gfp(:,j)'*Map_gfp(:,j)-(gamma(:,L(j))'*Map_gfp(:,j))^2);     
        end
        sigma_u = sigma_u/(frame*(num_channels-1)); 
    end
    if sigma_u < sigma
        sigma = sigma_u;
        Gamma = gamma;
    end
end

%% Segmentation Smoothing Algorithm
sigma_o = 0; sigma_u = 1000;
[~,L] = max((Map_gfp'*Gamma).^2,[],2);                 % compute L again        
Ka = L;

%% Calculate e
e = 0;
for j = 1:frame
    e = e + (Map_gfp(:,j)'*Map_gfp(:,j)-(Gamma(:,L(j))'*Map_gfp(:,j))^2);     
end
e = e/(frame*(num_channels-1));

times = 0;
while abs(sigma_o-sigma_u) > tol*sigma_u && times <  iter_thereshold
    sigma_o = sigma_u;
    
    % Compute Nbkt between a window whose length is 2b+1
    for t = 1+b:frame-b
        for k = 1:K
            Nbkt(t-b,k) = length(find(L(t-b:t+b) == k));
        end
        
        % argmin(k)
        S = zeros(K,1);
        for k = 1:K
            S(k) = (Map_gfp(:,t)'*Map_gfp(:,t)-(Gamma(:,k)'*Map_gfp(:,t))^2)/(2*e*(num_channels-1))-lamda*Nbkt(t-b,k);
        end
        [~,Min_index] = min(S);
        Ka(t) = Min_index; 
    end
    L = Ka;
    
    % Compute new sigma_u
    sigma_u = 0;
    for j = 1:frame
        sigma_u = sigma_u + (Map_gfp(:,j)'*Map_gfp(:,j)-(Gamma(:,L(j))'*Map_gfp(:,j))^2);     
    end
    sigma_u = sigma_u/(frame*(num_channels-1));
    times = times + 1;
end

%% Compute Final Params
Kappa = L;
A = 1:K;
a = Map_gfp' * Gamma;
for i =1:frame
    index =  A~=L(i);   
    a(i,index) = 0;
end

R_square = GEV; 
temp= 0;
for j = 1:frame
    temp = temp + (Map_gfp(:,j)'*Map_gfp(:,j)-(Gamma(:,L(j))'*Map_gfp(:,j))^2);     
end
R_cv = temp/(frame*(num_channels-1));
R_cv = R_cv*((num_channels-1)/(num_channels-1-K))^2;

%% Compute KL 
D_KL = zeros(K,1);
W = 0;
for i = 1:K
    index = find(Kappa == i);
    num_map(i) = length(index);
    D_KL(i) = 0;
    for j = 1: num_map(i)
        for k = j:num_map(i)   
            D_KL(i) = D_KL(i) + sum((Map_gfp(:,index(j))-Map_gfp(:,index(k))).^2);
        end    
    end
    W = W+D_KL(i)/(2*num_map(i));
end
M = W*K^(2/num_channels);

end
