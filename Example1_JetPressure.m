%% EXAMPLE 1: Compute CS-SPOD on the pressure of the St_f=1.5, a_0/U_j=10% forced turbulent jet.

% Clear variables
clc; clear;

% Make plotting pretty
set(0,'defaulttextinterpreter', 'latex');
set(0,'defaultlinelinewidth', 6 )
set(0,'defaultAxesFontSize', 40)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Add paths
addpath ./data/jet/
addpath ./src/

% Load in the data
dat     = load('./data/jet/JetDataForced.mat');
X       = dat.P;

% Define paramerters
dt      = dat.dt;                               % Timestep 
nt0     = dat.nt0;                              % Snapshots per oscillation
a0      = dat.a0;                               % Fundamental frequency
nt      = size(X, 1);                           % Total number of snapshots
Ma      = dat.Ma;                               % Mach number

% Load in the grid...
x       = dat.X;                                % X-Grid
r       = dat.R;                                % R-Grid

% Get the weighting matrix
intWeight       = trapzWeightsPolar(r(:,1),x(1,:));

% Define special estimation parameters
c       = 20;                                   % Scale factor on window length
LWindow = nt0*c;                                % Length of the window
WindowM = hamming(LWindow);                     % Load in the hamming window

% Get default parameters
opts    = defaultparams(X, nt0);                % Get the default parameters

% Define out own parameters
opts.Window         = WindowM;                  % Use our custom window
opts.Weight         = intWeight;                % Define our weighting matrix
opts.FracOvlp       = 0.67;                     % Fraction overlap
opts.singlef        = 0;                        % Compute the double sided spectrum
opts.nsave          = 2;                        % Save 5 dominant CS-SPOD modes
opts.blockwisemean  = false;                    % Do not remove the blockwise mean
opts.Nf             = 'max';                    % Use the max possible Nf (could also be a number >= 1)
opts.quiet          = 0;                        % Display the output

%% Compute CS-SPOD: Method 1

% Get the CS-SPOD modes
% L = [ngamma, nmodes];
% P = [ngamma, nx1, nx2, nx3, ..., nf_total, nmodes] where nx1, nx2 are the
% 1'st, 2'nd, etc spatial dimensions of the input matrix X.
% The flattened eigenvectors are = [ngamma, nx1*nx2*nx3*...*nf_total, nmodes]

% This is the efficient method as presented in the paper (algo 2).
% Memory is approx Mem(complex(X))*(1/(1 - Ovlp)+1)
[L, P, gamma, Nf, eTot] = csspod(X, a0, nt0, dt, opts);

%% Compute CS-SPOD: Method (Modified to reduce memory usage)
% This is the efficient method as presented in the paper (algo 2) but is slightly modified to
% reduce the memory usage. However, this method is much slower as you computer the fft nGamma times but the
% memory requirement is approx Mem(complex(X))*((Nf*2+1)/(nDFT*(1 - FracOvlp))+1) instead of
% Mem(complex(Q))*(1/(1 - FracOvlp)+1). Since (Nf*2+1)/(nDFT*(1 - FracOvlp))  << 1/(1 - FracOvlp)
% the memory requirement is much lower.
% [L, P, gamma, Nf, eTot] = csspod_memeff(X, a0, nt0, dt, opts);

%% Compute CS-SPOD: Method 2
% This is the original method (algo 3). This is not memory or computational
% time efficient and it is not reccomended to be used.
% [L, P, gamma, Nf, eTot] = csspod_original(X, a0, nt0, dt, opts);

%% Plot the eigenspectrum
figure(99)
clf
semilogy(gamma/Ma, L(:, 1), 'black-o'); hold on
semilogy(gamma/Ma, L(:, 2), 'blue-o');

axis square
xlabel('$\gamma$', 'interpreter', 'latex');
ylabel('$\lambda_j$', 'interpreter', 'latex');
xlim([-0.76 0.76])
legend('$\lambda_1$', '$\lambda_2$', 'interpreter', 'latex');
%% Plot the CS-SPOD modes (note that these scripts assume that the eigenvectors have been inflated)
iGam = 20; % Pick the gamma value to plot
iMode = 1; % Pick with mode to plot
disp(['Current gamma is : ', num2str(gamma(iGam)/Ma)]);

% Get the current mode
cMode = squeeze(P(iGam, :, :, :, iMode));
cGam = gamma(iGam);

% Define the time vector you are plotting
tM = 0:dt:(nt0*4-1)*dt; % Plot it for 10 oscillations
pM = tM/(nt0*dt)*2*pi; % Convert to phase
pM_indx =  mod(round(tM/dt), nt0) + 1;

% Get the time domain version of the dominant CS-SPOD mode
P_time = genTimeMode(cMode, cGam, a0, tM, Nf);

% Plot the real and absolute value of the mode for gamma = 0.05
figure(97)
disp('This video shows how the mode evolves in time')
for p = 1:length(tM)
    clf
    % Plot the current mode
    pcolor(x, r, real(P_time(:, :, p))); axis equal; shading interp; caxis([-1 1]*max(abs(P_time(:)))*0.5);
    axis([0 7 0 2]); colormap('jet');
    hold on
    
    % Overlay the phase-dependent mean velocity
    UxMeanCur = squeeze(dat.UxMean(pM_indx(p), :, :));
    contour(x, r, UxMeanCur, [1 1]*0.25, 'black-', 'lineWidth', 3);
    contour(x, r, UxMeanCur, [1 1]*0.75, 'black--', 'lineWidth', 3);
    xlabel('x');
    ylabel('r');
    
    % Pause before going to the next frame
    pause(0.1)
end

