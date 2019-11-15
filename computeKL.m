function KL = computeKL(M)
%% Compute Krzanowski-Lai criterion
%
% inputs:   1)M            : KL criterion params obtained from M-Kmeans
%
% ouputs:   1)KL           : Krzanowski-Lai criterion value
%  
% Anthor: Allard Shi, Xian Jiaotong University

num = size(M,2);
KL = zeros(1,num-2);
for i = 1:num-2
    if M(i) < M(i+1)
        KL(i) = 0;
    else
    KL(i) = 1+(M(i+2)-2*M(i+1))/M(i);
        if KL(i) < 0
            KL(i) = 0;
        end
    end
end
end
