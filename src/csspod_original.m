function [L, P, gammaM, Nf, eTot] = csspod_original(X, a0, nt0, dt, opts)
%%   CS-SPOD: Cyclostationary Spectral proper orthogonal decomposition
%   [L,P,gamma, Nf, eTot] = csspod_original(X, a0, nt0, dt, opts)
%   This algorithm efficiently computes CS-SPOD modes using algorithm 2
%   from [3]. This is an inefficient algorithm that should not be used 
%   unless absolutely required.CS-SPOD is an extension of SPOD [4] to 
%   flows with peirodic statistics. It is a multidimensional version of 
%   [3] using the method-of-snapshot technique for computational feasibility. 
%   CS-SPOD has been applied to a turbulent jet in [1, 2]. 


%   Inputs:
%   X (matrix): Data Matrix. First dimension is time and X can have any number
%               of spatial (or variable) dimensions
%   a0 (scalar): Fundamental frequency
%   nt0 (int): Number of snapshots per forcing oscillation
%   dt (double): Timestep between the snapshots. Note, dt = 1/(a0*nt0)
%   opts: A structure that contains parameters required to run CS-SPOD
%   opts.Window (vec): Windowing function to reduce spectral and cyclic aliasing
%   opts.Weight (vec/matrix): Weighting matrix
%   opts.FracOvlp (scalar betweem 0 - 1): Fraction that each block overlaps (to reduce variance
%   and spectral/cyclic leakage). Min, Max value = [0, 1]
%   opts.singlef (true, false): Single sided or double sided spectrum calculation
%   opts.nsave (int): How many dominant CS-SPOD modes to save
%   opts.blockwisemean (true, false): Subtract the blockwise mean
%   opts.Nf ('max' or scalar >= 1): Number of connected frequencies (systems) to use
%   opts.quiet (true, false): Display progress information
%   opts.inflate (true, false): Do you wish to inflate the final
%   eigenvectors, i.e. convert from a [nx1*nx2*nx3*...*nf_total, nsave]
%   matrix to [ngamma, nx1, nx2, nx3, ..., nf_total, nsave]?

%   Outputs:
%   L    : Eigenvalues ([ngamma, nsave] in size)
%   P    : Eigenvectors ([nx1*nx2*nx3*...*nf_total, nsave] or [ngamma, nx1, nx2, nx3, ..., nf_total, nsave])
%   gamma: Frequencies ([ngamma, 1])
%   Nf   : Number of connected systems (scalar)
%   eTot : Total energy per frequnecy in gamma ([ngamma, 1])
%
%
% % References:
% [1] Heidt, L. and Colonius, T., 2023. Spectral proper orthogonal decomposition of harmonically forced turbulent flows. arXiv preprint arXiv:2305.05628.
% [2] Heidt, L., Colonius, T., Nekkanti, A., Schmidt, O.T., Maia, I. and Jordan, P., 2023. Cyclostationary analysis of forced turbulent jets. In AIAA AVIATION 2023 Forum (p. 3652).
% [3] Kim, K.Y., North, G.R. and Huang, J., 1996. EOFs of one-dimensional cyclostationary time series: Computations, examples, and stochastic modeling. Journal of Atmospheric Sciences, 53(7), pp.1007-1017.
% [4] Lumley, J. L., Stochastic tools in turbulence, Academic Press, 1970

%% Please cite the following publications if you are using this code/method for your own research 
% Heidt, L. and Colonius, T., 2023. Spectral proper orthogonal decomposition of harmonically forced turbulent flows. arXiv preprint arXiv:2305.05628.
% @article{heidt2023spectral,
%   title={Spectral proper orthogonal decomposition of harmonically forced turbulent flows},
%   author={Heidt, Liam and Colonius, Tim},
%   journal={arXiv preprint arXiv:2305.05628},
%   year={2023}
% }
%
% Heidt, L., Colonius, T., Nekkanti, A., Schmidt, O.T., Maia, I. and Jordan, P., 2023. Cyclostationary analysis of forced turbulent jets. In AIAA AVIATION 2023 Forum (p. 3652).
% @inproceedings{heidt2023cyclostationary,
%   title={Cyclostationary analysis of forced turbulent jets},
%   author={Heidt, Liam and Colonius, Tim and Nekkanti, Akhil and Schmidt, Oliver T and Maia, Igor and Jordan, Peter},
%   booktitle={AIAA AVIATION 2023 Forum},
%   pages={3652},
%   year={2023}
% }

% Authors: Liam Heidt (lheidt@caltech.edu), Tim Colonius
% Affiliation: California Institute of Technology

%% License
% MIT License
% Copyright (c) 2023 Liam Heidt, Tim Colonius (California Institute of Technology)

%% Last revision: 10-Nov-2023 (Liam Heidt)

%% Unpack the structure to improve readability
Window                  = opts.Window;
Weight                  = opts.Weight;
FracOvlp                = opts.FracOvlp;
singlef                 = opts.singlef;
nsave                   = opts.nsave;
blockwisemean           = opts.blockwisemean;
Nf                      = opts.Nf;
quiet                   = opts.quiet;
delf                    = 1/(length(Window)*dt);

tic
if ~quiet
    disp('********************************************************************');
    disp('Starting CS-SPOD computation using original (inefficient) algorithm');
    disp('********************************************************************');
end

% Output the parameters
if ~quiet
    disp('    - Displaying parameters');
    disp(['        - Window length: ', num2str(length(Window)), ' (c = ', num2str(length(Window)/nt0), ', \Delta gamma = ', num2str(delf, 2), ')']);
    disp(['        - Fraction Overlap: ', num2str(FracOvlp)]);
    if singlef
        disp('        - Single sided spectrum computation');
    else
        disp('        - Double sided spectrum computation');
    end
    disp(['        - Computing ', num2str(nsave), ' modes']);
    if blockwisemean
        disp('        - Subtracting blockwise mean');
    else
        disp('        - Not subtracting blockwise mean');
    end
    if strcmp(Nf, 'max')
        disp('        - Using max number of linked frequnecies')
    else
        disp(['        - Using Nf = ', num2str(Nf), ' linked frequencies'])
    end
end
%% Flatten all input matrices (data and weight)
if ~quiet; disp(['    - Flattening input matrices']); end

% First dimension is always time
nt = size(X, 1);
% 2nd to end dimensions are spatial
XSize = size(X); XSize = XSize(2:end);
WSize = size(Weight); WSize = WSize(1:length(XSize));
if ~isequal(XSize, WSize)
    error('Spatial size of data matrix (X) and weights (W) is not consistent');
end

% Total spatial dimension
nx                      = prod(XSize);

% Flatten the matrices
X       = reshape(X, [nt, nx]);
Weight  = reshape(Weight, [nx, 1]);

%% Set up the frequency axis
if ~quiet; disp('    - Setting up frequency index'); end

% Length of the window
nDFT = length(Window);
if nDFT < 4
    error('Spectral estimation parameters not meaningful: nDFT is too short.');
end

% Compute how many connected frequnecies to use
if strcmp(Nf, 'max')
    max_f = 1/dt/2;
    Nf = floor((max_f - a0/2)/a0);
end
% Total number of connected frequencies
nf_total = 2*Nf+1;

% Frequency shifts of the data (as a multiple of alpha)
nfM = -Nf:1:Nf;

% obtain frequency axis (equation 3.16)
f = (0:nDFT-1)/(nDFT*dt);
if mod(nDFT,2)==0
    f(nDFT/2+1:end) = f(nDFT/2+1:end)-1/dt;
else
    f((nDFT+1)/2+1:end) = f((nDFT+1)/2+1:end)-1/dt;
end


% Compute the index of the values of gamma (equation 3.16)
gammaM_k = 1:floor(a0*nDFT*dt/2)+1; % k <= floor(c/2)+1
gammaM = (gammaM_k-1)/(nDFT*dt);

if ~singlef
    gammaM_k2 = nDFT-ceil(a0*nDFT*dt/2)+1+1:nDFT; % N_f - ceil(c/2)+1 < k <= Nf
    gammaM_k = [gammaM_k, gammaM_k2];
    gammaM   = [gammaM, (gammaM_k2 - 1 - nDFT)/(nDFT*dt)];
end

% Get the gamma values from the frequency vector and ensure it matches the analytically expected ones
gammaM_V2 = f(gammaM_k);
numGamma = length(gammaM);
if max(gammaM_V2 - gammaM) > 1e-10
    error('True gamma values and expected gamma values do not match');
end

%% Compute the frequency data matrix
if ~quiet; disp('    - Computing the data matrix'); end

% Compute the number of snapshots to overlap by
nOvlp = ceil(nDFT*FracOvlp); % Rounding up to ensure we have a minimum overlap

% Compute the number of realizations we will have in total
nBlks           = floor((nt-nOvlp)/(nDFT-nOvlp));
if nBlks < 3
    error('Spectral estimation parameters not meaningful: nBlks is too small.');
end

% Init the data matrix with the same type as the orignal matrix
dataM = zeros(nDFT, nx, nBlks, class(X));

% Get the time vector
tM = (0:1:nt)*dt;

% Init the time matrix
tM_M = zeros(nDFT, nBlks);

% Loop through all of the blocks
for iBlk = 1:nBlks
    
    % Calculate the offset required for the current block (this is the
    % offset required to extract the data for the current block as well as
    % the phase-shift required for the efificent method)
    offset = (iBlk-1)*(nDFT-nOvlp);
    
    % Get the index of the current data matrix
    timeIndx = (1:nDFT) + offset;
    
    % Get the current data
    curDat = X(timeIndx, :);
    curTime = tM(timeIndx);
    
    % Remove the block mean if required
    if blockwisemean
        curDat = curDat - mean(curDat, 1);
    end
    
    % Save the data to the data matrix
    dataM(:, :, iBlk) = curDat;
    tM_M(:, iBlk) = curTime;
end

%% Now compute the frequency data matrix
if ~quiet; disp('    - Computing the frequency data matrix'); end

% Initialize the matrices that will store the final eigenvalues and
% eigenvectors
numGam  = length(gammaM);

% Init the data matrix with the same type as the orignal matrix
dataM_Shifted = zeros(numGam, nx, nf_total, nBlks, class(X)); dataM_Shifted = complex(dataM_Shifted, 0);

% Loop over all gamma values
for iNf = 1:nf_total
    
    if ~quiet; disp(['        - Cyclic frequency: ', num2str(iNf), ' out of ', num2str(nf_total)]); end
    
    % Get the current alpha
    alpha = nfM(iNf)*a0;
    
    % Loop over all blocks
    for iBlk = 1:nBlks
        % Freq shift the data and then take the fourier transform
        cDat = fft((dataM(:, :, iBlk).*exp(-1i*2*pi*alpha*tM_M(:, iBlk))).*Window, [], 1);
        dataM_Shifted(:, :, iNf, iBlk) = cDat(gammaM_k, :);
    end
    
end

% Consturct the concatenated frequency-data matrix
dataM_Shifted = reshape(dataM_Shifted, [numGam, nx*nf_total, nBlks]);

%% Loop over all the values of gamma and compute CS-SPOD for each value of gamma
if ~quiet; disp('    - Computing CS-SPOD for each gamma'); end

% Initialize the matrices that will store the final eigenvalues and
% eigenvectors
P       = zeros(numGam, nx*nf_total, nsave, class(X)); P = complex(P, 0);
L       = zeros(numGam, nsave, class(X)); L = complex(L, 0);
eTot    = zeros(numGam, 1, class(X));

for iGam = 1:numGam
    if ~quiet; disp(['        - Gamma index: ', num2str(iGam), ' out of ', num2str(numGamma)]); end
    
    
    %% Compute CS-SPOD using the method of snapshot technique
    
    % Generate the concatenated frequency-data matrix (equation 3.14b)
    % using equation 3.23b to efficiently compute it
    Qtilde = squeeze(dataM_Shifted(iGam, :, :));
    
    % Determine kappa (normalization constant)
    kappa = dt/(nBlks*norm(Window)^2);
    
    % Scale the concatenated frequency-data matrix
    Qtilde = Qtilde*sqrt(kappa);
    
    % Determine the total weight matrix (repeat Weight along the diagonal
    % nf_total times)
    WeightT = repmat(Weight, nf_total, 1);                            
    
    % Compute M  = Qtilde'*W*Qtilde (equation 3.20)
    M = double(Qtilde'*(WeightT.*Qtilde));
    
    % Compute the eigendecomposition
    [Theta,Lambda]      = eigs(M,  nsave, 'largestabs', 'Tolerance', 1e-12);
    Lambda = diag(Lambda);
    Psitilde = Qtilde*Theta*diag(1./sqrt(Lambda));
    
    % Save the eigenvalues and vectors to L and P
    L(iGam, :) = real(Lambda);
    P(iGam, :, :) = Psitilde;
    
    % Other output...
    eTot(iGam) = sum(abs(diag(M)));
    
end

%% Shift the frequency to go from lowest to highest
if ~singlef
    pos_f    = 1:floor(numGamma/2)+1;
    neg_f    = floor(numGamma/2)+2:numGamma;
    gammaM   = gammaM([neg_f, pos_f]);
    L        = L([neg_f, pos_f], :);
    P        = P([neg_f, pos_f], :);
    eTot     = eTot([neg_f, pos_f]);
end


%% Inflate the matrices
if opts.inflate
    if ~quiet; disp('    - Inflating CS-SPOD output modes'); end
    P = reshape(P, [numGam, XSize ,nf_total, nsave]);
end

if ~quiet
    disp('*******************************************************************');
    disp(['          Finished CS-SPOD computation in ', num2str(toc), ' seconds          ']);
    disp('*******************************************************************');
end
end
