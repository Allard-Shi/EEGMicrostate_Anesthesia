function C = compute_spatial_correlation(tseries,ms)
%% Compute spatial correlation
%
% inputs:   1)tseries  : Time seies map (M*N matrix)
%                        M - number of electrodes  N - number of time frame
%           2)ms       : Microstate     (M*N matrix)
%                        M - number of electrodes  N - number of microstate
%
% ouputs:   1)C        : Correlation    (M*N matrix)
%                        M - number of time frame  N - number of microstate
%  

%% Spatial correlation
[~,frame] = size(tseries);
[~,num_ms] = size(ms);
C = zeros(frame,num_ms);
for i = 1:num_ms
    for j = 1:frame
        C(j,i) = sum(tseries(:,j).*ms(:,i))/sqrt(sum(tseries(:,j).^2)*sum(ms(:,i).^2));
		
		% absolute correlation 
        C(j,i) = abs(C(j,i));
    end
end

end