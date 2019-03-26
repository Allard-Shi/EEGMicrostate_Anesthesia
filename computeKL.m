function KL = computeKL(W,num_electrode)
%% Compute Krzanowski-Lai criterion
%
% inputs:   1)W            : KL criterion params (Dispersion)
%
% ouputs:   1)KL           : Krzanowski-Lai criterion value
%           2)num_electrode: Number of electrodes
%  
% Anthor: Allard Shi, Xian Jiaotong University

num = size(W,1);
KL = zeros(1,num-2);
for i =1:num
    M = W*i^(2/num_electrode);
end
for i = 1:num-2
    KL(i) = 1+(M(i+2)-2*M(i+1))/M(i);
end
end
