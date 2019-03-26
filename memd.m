function q = memd(x, varargin)
%
%
% function MEMD applies the "Multivariate Empirical Mode Decomposition" algorithm (Rehman and Mandic, Proc. Roy. Soc A, 2010)
% to multivariate inputs. We have verified this code by simulations for signals containing 3-16 channels.
%
% Syntax:
%
% imf = MEMD(X)
%   returns a 3D matrix 'imf(N,M,L)' containing M multivariate IMFs, one IMF per column, computed by applying
%   the multivariate EMD algorithm on the N-variate signal (time-series) X of length L.
%    - For instance, imf_k = IMF(k,:,:) returns the k-th component (1 <= k <= N) for all of the N-variate IMFs.
%
%   For example,  for hexavariate inputs (N=6), we obtain a 3D matrix IMF(6, M, L)
%   where M is the number of IMFs extracted, and L is the data length.
%
% imf = MEMD(X,num_directions)
%   where integer variable num_directions (>= 1) specifies the total number of projections of the signal
%     - As a rule of thumb, the minimum value of num_directions should be twice the number of data channels,
%     - for instance, num_directions = 6  for a 3-variate signal and num_directions= 16 for an 8-variate signal
%   The default number of directions is chosen to be 64 - to extract meaningful IMFs, the number of directions
%   should be considerably greater than the dimensionality of the signals
%
% imf = MEMD(X,num_directions,'stopping criteria')
%   uses the optional parameter 'stopping criteria' to control the sifting process.
%    The available options are
%      -  'stop' which uses the standard stopping criterion specified in [2]
%      -  'fix_h' which uses the modified version of the stopping criteria specified in [3]
%    The default value for the 'stopping criteria' is 'stop'.
%
%  The settings  num_directions=64 and 'stopping criteria' = 'stop' are defaults.
%     Thus imf = MEMD(X) = MEMD(X,64) = MEMD(X,64,'stop') = MEMD(X,[],'stop'),
%
% imf = MEMD(X, num_directions, 'stop', stop_vec)
%   computes the IMFs based on the standard stopping criterion whose parameters are given in the 'stop_vec'
%     - stop_vec has three elements specifying the threshold and tolerance values used, see [2].
%     - the default value for the stopping vector is   step_vec = [0.075 0.75 0.075].
%     - the option 'stop_vec' is only valid if the parameter 'stopping criteria' is set to 'stop'.
%
% imf = MEMD(X, num_directions, 'fix_h', n_iter)
%   computes the IMFs with n_iter (integer variable) specifying the number of consecutive iterations when
%   the number of extrema and the number of zero crossings differ at most by one [3].
%     - the default value for the parameter n_iter is set to  n_iter = 2.
%     - the option n_iter is only valid if the parameter  'stopping criteria' = 'fix_h'
%
%
% This code allows to process multivaraite signals having 3-16 channels, using the multivariate EMD algorithm [1].
%   - to perform EMD on more than 16 channels, modify the variable 'Max_channels' on line 510 in the code accordingly.
%   - to process 1- and 2-dimensional (univariate and bivariate) data using EMD, we recommend the toolbox from
%                 http://perso.ens-lyon.fr/patrick.flandrin/emd.html
%
% Acknowledgment: Part of this code is based on the bivariate EMD code, publicly available from
%                 http://perso.ens-lyon.fr/patrick.flandrin/emd.html. We would also like to thank 
%                 Anh Huy Phan from RIKEN for helping us in optimizing the code and making it computationally efficient. 
%
%
% Copyright: Naveed ur Rehman and Danilo P. Mandic, Oct-2009
%
%
% [1]  Rehman and D. P. Mandic, "Multivariate Empirical Mode Decomposition", Proceedings of the Royal Society A, 2010
% [2]  G. Rilling, P. Flandrin and P. Goncalves, "On Empirical Mode Decomposition and its Algorithms", Proc of the IEEE-EURASIP
%      Workshop on Nonlinear Signal and Image Processing, NSIP-03, Grado (I), June 2003
% [3]  N. E. Huang et al., "A confidence limit for the Empirical Mode Decomposition and Hilbert spectral analysis",
%      Proceedings of the Royal Society A, Vol. 459, pp. 2317-2345, 2003

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Usage %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case 1:

% inp = randn(1000,3);
% imf = memd(inp);
% imf_x = reshape(imf(1,:,:),size(imf,2),length(inp)); % imfs corresponding to 1st component
% imf_y = reshape(imf(2,:,:),size(imf,2),length(inp)); % imfs corresponding to 2nd component
% imf_z = reshape(imf(3,:,:),size(imf,2),length(inp)); % imfs corresponding to 3rd component


% Case 2:

% load syn_hex_inp.mat
% imf = memd(s6,256,'stop',[0.05 0.5 0.05])

global N N_dim;
[x, seq, t, ndir, N_dim, N, sd, sd2, tol, nbit, MAXITERATIONS, stop_crit, stp_cnt] = set_value(x, nargin, varargin{:});

r=x; n_imf=1;q = zeros(N_dim,1,N);
while ~stop_emd(r, seq, ndir)
    % current mode
    m = r;
     
    % computation of mean and stopping criterion
    if(strcmp(stop_crit,'stop'))
        [stop_sift,env_mean] = stop_sifting(m,t,sd,sd2,tol,seq,ndir);
    else
        counter=0;
        [stop_sift,env_mean,counter] = stop_sifting_fix(m,t,seq,ndir,stp_cnt,counter);
    end
    
    % In case the current mode is so small that machine precision can cause
    % spurious extrema to appear
    if (max(abs(m))) < (1e-10)*(max(abs(x)))
        if ~stop_sift
            warning('emd:warning','forced stop of EMD : too small amplitude')
        else
            disp('forced stop of EMD : too small amplitude')
        end
        break
    end
    
    % sifting loop
    while ~stop_sift && nbit<MAXITERATIONS
        %sifting
        m = m - env_mean;
        % computation of mean and stopping criterion
        if(strcmp(stop_crit,'stop'))
            [stop_sift,env_mean] = stop_sifting(m,t,sd,sd2,tol,seq,ndir);
        else
            [stop_sift,env_mean,counter] = stop_sifting_fix(m,t,seq,ndir,stp_cnt,counter);
        end
    
        nbit=nbit+1;
        
        if(nbit==(MAXITERATIONS-1) &&  nbit > 100)
            warning('emd:warning','forced stop of sifting : too many iterations');
        end
    end
    
    q(:,n_imf,:)=m';
    
    n_imf = n_imf+1;
    r = r - m;
    nbit = 0;
end
% Stores the residue
q(:,n_imf,:)=r';

%sprintf('Elapsed time: %f\n',toc);
end

%---------------------------------------------------------------------------------------------------
function stp = stop_emd(r, seq, ndir)
global N_dim;
ner = zeros(ndir,1);
for it=1:ndir
    if (N_dim~=3) % Multivariate signal (for N_dim ~=3) with hammersley sequence
        % Linear normalisation of hammersley sequence in the range of -1.00 - 1.00
        b=2*seq(1:end,it)-1;
        
        % Find angles corresponding to the normalised sequence
        tht = atan2(sqrt(flipud(cumsum(b(N_dim:-1:2).^2))),b(1:N_dim-1)).';
        % Find coordinates of unit direction vectors on n-sphere
        dir_vec(1:N_dim) = [1 cumprod(sin(tht))];
        dir_vec(1:N_dim-1) =  cos(tht) .*dir_vec(1:N_dim-1);

    else % Trivariate signal with hammersley sequence
        % Linear normalisation of hammersley sequence in the range of -1.0 - 1.0
        tt = 2*seq(1,it)-1;
        tt((tt>1))=1;
        tt((tt<-1))=-1;
        
        % Normalize angle from 0 - 2*pi
        phirad = seq(2,it)*2*pi;
        st = sqrt(1.0-tt*tt);
        
        dir_vec(1)=st * cos(phirad);
        dir_vec(2)=st * sin(phirad);
        dir_vec(3)=tt;
    end
    % Projection of input signal on nth (out of total ndir) direction
    % vectors
    y = r * dir_vec';
    % Calculates the extrema of the projected signal
    [indmin, indmax] = local_peaks(y);
    
    ner(it) = length(indmin) + length(indmax);
end

% Stops if the all projected signals have less than 3 extrema
stp = all(ner < 3);
end

%---------------------------------------------------------------------------------------------------
% computes the mean of the envelopes and the mode amplitude estimate
function [env_mean,nem,nzm,amp] = envelope_mean(m,t,seq,ndir) %new
global N N_dim;
NBSYM = 2;
count=0;

env_mean=zeros(length(t),N_dim);
amp = zeros(length(t),1);
nem = zeros(ndir,1);nzm = zeros(ndir,1);
for it=1:ndir
    if (N_dim ~=3) % Multivariate signal (for N_dim ~=3) with hammersley sequence
        % Linear normalisation of hammersley sequence in the range of -1.00 - 1.00
        b=2*seq(1:end,it)-1;
        % Find angles corresponding to the normalised sequence
        tht = atan2(sqrt(flipud(cumsum(b(N_dim:-1:2).^2))),b(1:N_dim-1)).';
        % Find coordinates of unit direction vectors on n-sphere
        dir_vec(1:N_dim) = [1 cumprod(sin(tht))];
        dir_vec(1:N_dim-1) =  cos(tht) .*dir_vec(1:N_dim-1);
    else % Trivariate signal with hammersley sequence
        % Linear normalisation of hammersley sequence in the range of -1.0 - 1.0
        tt = 2*seq(1,it)-1;
        tt((tt>1))=1;
        tt((tt<-1))=-1;
        
        % Normalize angle from 0 - 2*pi
        phirad = seq(2,it)*2*pi;
        st = sqrt(1.0-tt*tt);
        
        dir_vec(1)=st * cos(phirad);
        dir_vec(2)=st * sin(phirad);
        dir_vec(3)=tt;
    end
    
    % Projection of input signal on nth (out of total ndir) direction vectors
    y(1:N)  = dir_vec * m(1:N,:)';
    
    % Calculates the extrema of the projected signal
    [indmin, indmax] = local_peaks(y);
    
    
    nem(it) = length(indmin) + length(indmax);
    
    indzer = zero_crossings(y);
    nzm(it) = length(indzer);
    
    [tmin,tmax,zmin,zmax,mode] = boundary_conditions(indmin,indmax,t,y,m,NBSYM);
    
    % Calculate multidimensional envelopes using spline interpolation
    % Only done if number of extrema of the projected signal exceed 3
    if(mode)
        env_min = spline(tmin,zmin.',t).';
        env_max = spline(tmax,zmax.',t).';
        
        amp = amp + sqrt(sum((env_max-env_min).^2,2))/2;
        env_mean = env_mean + (env_max+env_min)/2;
    else % if the projected signal has inadequate extrema
        count=count+1;
    end
end
if(ndir>count)
    env_mean = env_mean/(ndir-count);
    amp = amp/(ndir-count);
else
    env_mean = zeros(N,N_dim);
    amp = zeros(N,1);
    nem = zeros(1,ndir);
end
end

%-------------------------------------------------------------------------------
% Stopping criterion
function [stp,env_mean] = stop_sifting(m,t,sd,sd2,tol,seq,ndir)
global N N_dim;
try
    [env_mean,nem,nzm,amp] = envelope_mean(m,t,seq,ndir);
    sx = sqrt(sum(env_mean.^2,2));
    if(amp) % something is wrong here
        sx = sx./amp;
    end
    stp = ~((mean(sx > sd) > tol | any(sx > sd2)) & any(nem > 2));
catch
    env_mean = zeros(N,N_dim);
    stp = 1;
end
end

function [stp,env_mean,counter]= stop_sifting_fix(m,t,seq,ndir,stp_count,counter)
global N N_dim;
try
    [env_mean,nem,nzm] = envelope_mean(m,t,seq,ndir);
    if (all(abs(nzm-nem)>1))
        stp = 0;
        counter = 0;
    else
        counter = counter+1;
        stp = (counter >= stp_count);
    end
catch
    env_mean = zeros(N,N_dim);
    stp = 1;
end
end

%---------------------------------------------------------------------------------------
% defines new extrema points to extend the interpolations at the edges of the
% signal (mainly mirror symmetry)
function [tmin,tmax,zmin,zmax,mode] = boundary_conditions(indmin,indmax,t,x,z,nbsym)

lx = length(x);
if (length(indmin) + length(indmax) < 3)
    mode = 0;
    tmin=NaN;tmax=NaN;zmin=NaN;zmax=NaN;
    return
else
    mode=1; %the projected signal has inadequate extrema
end
% boundary conditions for interpolations :
if indmax(1) < indmin(1)
    if x(1) > x(indmin(1))
        lmax = fliplr(indmax(2:min(end,nbsym+1)));
        lmin = fliplr(indmin(1:min(end,nbsym)));
        lsym = indmax(1);
    else
        lmax = fliplr(indmax(1:min(end,nbsym)));
        lmin = [fliplr(indmin(1:min(end,nbsym-1))),1];
        lsym = 1;
    end
else
    
    if x(1) < x(indmax(1))
        lmax = fliplr(indmax(1:min(end,nbsym)));
        lmin = fliplr(indmin(2:min(end,nbsym+1)));
        lsym = indmin(1);
    else
        lmax = [fliplr(indmax(1:min(end,nbsym-1))),1];
        lmin = fliplr(indmin(1:min(end,nbsym)));
        lsym = 1;
    end
end

if indmax(end) < indmin(end)
    if x(end) < x(indmax(end))
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
        rmin = fliplr(indmin(max(end-nbsym,1):end-1));
        rsym = indmin(end);
    else
        rmax = [lx,fliplr(indmax(max(end-nbsym+2,1):end))];
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
        rsym = lx;
    end
else
    if x(end) > x(indmin(end))
        rmax = fliplr(indmax(max(end-nbsym,1):end-1));
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
        rsym = indmax(end);
    else
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
        rmin = [lx,fliplr(indmin(max(end-nbsym+2,1):end))];
        rsym = lx;
    end
end
tlmin = 2*t(lsym)-t(lmin);
tlmax = 2*t(lsym)-t(lmax);
trmin = 2*t(rsym)-t(rmin);
trmax = 2*t(rsym)-t(rmax);

% in case symmetrized parts do not extend enough
if tlmin(1) > t(1) || tlmax(1) > t(1)
    if lsym == indmax(1)
        lmax = fliplr(indmax(1:min(end,nbsym)));
    else
        lmin = fliplr(indmin(1:min(end,nbsym)));
    end
    if lsym == 1
        error('bug')
    end
    lsym = 1;
    tlmin = 2*t(lsym)-t(lmin);
    tlmax = 2*t(lsym)-t(lmax);
end

if trmin(end) < t(lx) || trmax(end) < t(lx)
    if rsym == indmax(end)
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
    else
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
    end
    if rsym == lx
        error('bug')
    end
    rsym = lx;
    trmin = 2*t(rsym)-t(rmin);
    trmax = 2*t(rsym)-t(rmax);
end
zlmax =z(lmax,:);
zlmin =z(lmin,:);
zrmax =z(rmax,:);
zrmin =z(rmin,:);

tmin = [tlmin t(indmin) trmin];
tmax = [tlmax t(indmax) trmax];
zmin = [zlmin; z(indmin,:); zrmin];
zmax = [zlmax; z(indmax,:); zrmax];
end

function [indmin, indmax] = local_peaks(x)
if(all(x < 1e-6))
    x=zeros(1,length(x));
end
m = length(x);
% Calculates the extrema of the projected signal
% Difference between subsequent elements:
dy = diff(x); a = find(dy~=0);
lm = find(diff(a)~=1) + 1;
d = a(lm) - a(lm-1);
a(lm) = a(lm) - floor(d/2);
a(end+1) = m;
ya  = x(a);

if(length(ya) > 1)
    % Maxima
    [pks_max,loc_max]=peaks(ya);
    % Minima
    [pks_min,loc_min]=peaks(-ya);
    
    if(~isempty(pks_min))
        indmin = a(loc_min);
    else
        indmin = NaN;
    end
    if(~isempty(pks_max))
        indmax = a(loc_max);
    else
        indmax = NaN;
    end
else
    indmin=NaN;
    indmax=NaN;
end
end

function [pks_max,locs_max] =peaks(X)
dX = sign(diff(X));
locs_max = find((dX(1:end-1) >0) &  (dX(2:end) <0)) + 1;
pks_max = X(locs_max);
end


function indzer = zero_crossings(x)
indzer = find(x(1:end-1).*x(2:end)<0);
if any(x == 0)
    iz = find( x==0 );
    if any(diff(iz)==1)
        zer = x == 0;
        dz = diff([0 zer 0]);
        debz = find(dz == 1);
        finz = find(dz == -1)-1;
        indz = round((debz+finz)/2);
    else
        indz = iz;
    end
    indzer = sort([indzer indz]);
end
end

function seq = hamm(n,base)
seq = zeros(1,n);
if ( 1 < base )
    seed = 1:1:n;
    base_inv = inv(base);
    while ( any ( seed ~= 0 ) )
        digit = mod (seed(1:n), base);
        seq = seq + digit * base_inv;
        base_inv = base_inv / base;
        seed = floor (seed / base );
    end
else
    temp = 1:1:n;
    seq = (mod(temp,(-base + 1 ))+0.5)/(-base);
end
end

function [q, seq, t, ndir, N_dim, N, sd, sd2, tol, nbit, MAXITERATIONS, stp_crit, stp_cnt] = set_value(q, narg, varargin)

error(nargchk(1,4,narg));
ndir = [];
stp_crit = [];
stp_vec = [];
stp_cnt  = [];
MAXITERATIONS  = [];
sd=[];
sd2=[];
tol=[];
prm= [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,...
      149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,...       
      307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,...                                                                                                             
      467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,...
      653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,...
      853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,...
      1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,... 
      1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,...
      1373,1381,1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,1523,1531,... 
      1543,1549,1553,1559,1567,1571,1579,1583,1597,1601,1607,1609,1613];
      
% Changes the input vector to double vector
q = double(q);

% Specifies maximum number of channels that can be processed by the code
% Its maximum possible value is 128.
Max_channels = 91;

if(narg==2)
    ndir=varargin{1};
end

if(narg==3)
    if(~isempty(varargin{1}))
        ndir=varargin{1};
    else
        ndir=182;
    end
    stp_crit=varargin{2};
end

if(narg==4 && strcmp(varargin{2},'fix_h'))
    if(isempty(varargin{1}))
        ndir=182;
        stp_crit=varargin{2};
        stp_cnt  = varargin{3};
    else
        ndir=varargin{1};
        stp_crit=varargin{2};
        stp_cnt  = varargin{3};
    end
elseif (narg==4 && strcmp(varargin{2},'stop'))
    if(isempty(varargin{1}))
        ndir=182;
        stp_crit=varargin{2};
        stp_vec=varargin{3};
    else
        ndir=varargin{1};
        stp_crit=varargin{2};
        stp_vec=varargin{3};
    end
elseif (narg==4 && ~xor(strcmp(varargin{2},'fix_h'),strcmp(varargin{2},'stop')))
    Nmsgid = generatemsgid('invalid stop_criteria');
    error(Nmsgid,'stop_criteria should be either fix_h or stop');
end

%%%%%%%%%%%%%% Rescale input signal if required
if (any(size(q)) == 0)
    datamsgid = generatemsgid('emptyDataSet');
    error(datamsgid,'Data set cannot be empty.');
end
if size(q,1) < size(q,2)
    q=q';
end

%%%%%%%%%%%% Dimension of input signal
N_dim = size(q,2);
if(N_dim < 3 || N_dim > Max_channels)
    error('Function only processes the signal having 3 and 16 channels.');
end

%%%%%%%%%%%% Length of input signal
N = size(q,1);

%%%%%%%%%%%%% Check validity of Input parameters
if ~isempty(ndir) && (~isnumeric(ndir) || ~isscalar(ndir) || any(rem(ndir,1)) || (ndir < 4))
    Nmsgid = generatemsgid('invalid num_dir');
    error(Nmsgid,'num_dir should be an integer greater than or equal to 6.');
end

if ~isempty(stp_crit) && (~ischar(stp_crit) || ~xor(strcmp(stp_crit,'fix_h'),strcmp(stp_crit,'stop')))
    Nmsgid = generatemsgid('invalid stop_criteria');
    error(Nmsgid,'stop_criteria should be either fix_h or stop');
end

if ~isempty(stp_vec) && (~isnumeric(stp_vec) || length(stp_vec)~=3 || ~strcmp(stp_crit,'stop'))
    Nmsgid = generatemsgid('invalid stop_vector');
    error(Nmsgid,'stop_vector should be an array with three elements e.g. default is [0.075 0.75 0.075] ');
end

if ~isempty(stp_cnt) && (~isnumeric(stp_cnt) || ~isscalar(stp_cnt) || any(rem(stp_cnt,1)) || (stp_cnt < 0) || ~strcmp(stp_crit,'fix_h'))
    Nmsgid = generatemsgid('invalid stop_count');
    error(Nmsgid,'stop_count should be a nonnegative integer');
end

if (isempty(ndir))
    ndir = 182; % default
end

if (isempty(stp_crit))
    stp_crit='stop'; % default
end

if (isempty(stp_vec))
    stp_vec=[0.075,0.75,0.075]; % default
end

if (isempty(stp_cnt))
    stp_cnt=2; % default
end

if(strcmp(stp_crit,'stop'))
    sd = stp_vec(1);
    sd2 = stp_vec(2);
    tol = stp_vec(3);
end

%%%%%%%%%%%%% Initializations for Hammersley function
base(1) = -ndir;

%%%%%%%%%%%%%% Find the pointset for the given input signal
if(N_dim==3)
    base(2) = 2;
    for it=1:N_dim-1
        seq(it,:) = hamm(ndir,base(it));
    end
  
else
    for iter = 2 : N_dim
        base(iter) = prm(iter-1);
    end
    
    for it=1:N_dim
        seq(it,:) = hamm(ndir,base(it));
    end
end

%%%%%%%%%%%% Define t
t=1:N;

% Counter
nbit=0;
MAXITERATIONS=1000; % default

% tic
end