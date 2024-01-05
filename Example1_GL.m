%% EXAMPLE 1: Compute CS-SPOD on Ginzburg Landau data.

% Clear variables
clc; clear;

% Make plotting nice
set(0,'defaulttextinterpreter', 'latex');
set(0,'defaultlinelinewidth', 6 )
set(0,'defaultAxesFontSize', 40)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

% Add path
addpath ./src/

% Load in the data
if isfile('./data/gl/GL_CS_Amu_0p4.mat') 
    dat     = load('./data/gl/GL_CS_Amu_0p4.mat');
else
    error('To run the examples, download data from: https://caltech.box.com/v/CSSPODGithubData');
end
X       = dat.solSaveM;                         % Data vector (nt, nx1, nx2, nx3)
a0      = dat.consts.omega;                     % Fundamental cyclic frequency (forcing freq)
nt      = size(X, 1);                           % Total number of snapshots
nt0     = dat.consts.np;                        % Number of snapshots per oscillation
dt      = 1/(a0*nt0);                           % Timestep
x       = dat.xgrid;                            % Grid

% Define special estimation parameters
c = 20;                                         % Scale factor on window length
LWindow = nt0*c+7;                                % Length of the window
WindowM = hamming(LWindow);                     % Load in the window

% Get default parameters
opts = defaultparams(X, nt0);

% Compute the actual weight matrix to use
dx = diff(x);
dx = [dx(1); dx; dx(end)];
Weight = 0.5*(dx(1:end-1) + dx(2:end));         % Matrix of (nx1, nx2, nx3)

% Define out own parameters
opts.Window         = WindowM;                  % Use our custom window
opts.Weight         = Weight;                   % Define our weighting matrix
opts.FracOvlp       = 0.75;                     % Frac overlap
opts.singlef        = 0;                        % Compute the double sided spectrum
opts.nsave          = 3;                        % Save 5 dominant CS-SPOD modes
opts.blockwisemean  = false;                    % Do not remove the blockwise mean
opts.Nf             = 'max';                    % Use the max possible Nf (could also be a number >= 1)
opts.quiet          = 0;                        % Display the output

%% Compute CS-SPOD: Method 1

% Get the CS-SPOD modes
% L = [ngamma, nmodes];
% P = [ngamma, nx1, nx2, nx3, ..., nf_total, nmodes] where nx1, nx2 are the
% 1'st, 2'nd, etc spatial dimensions of the input matrix X. 
% The flattened eigenvectors are = [ngamma, nx1*nx2*nx3*...*nf_total, nmodes]

%% Compute CS-SPOD: Method 1 (Modified to reduce memory usage)
% This is the efficient method as presented in the paper (algo 2). 
% Memory is approx Mem(complex(X))*(1/(1 - Ovlp)+1)
[L, P, gamma, Nf, eTot] = csspod(X, a0, nt0, dt, opts);

%% Compute CS-SPOD: Method 2
% This is the efficient method as presented in the paper (algo 2). 
% However, this method is slower as you computer the fft nGamma times but the
% memory requirement is approx Mem(complex(X))*((Nf*2+1)/(nDFT*(1 - FracOvlp))+1) instead of 
% Mem(complex(Q))*(1/(1 - FracOvlp)+1). Since (Nf*2+1)/(nDFT*(1 - FracOvlp))  << 1/(1 - FracOvlp)
% the memory requirement is much lower. 
% [L, P, gamma, Nf, eTot] = csspod_memeff(X, a0, nt0, dt, opts); 


% This is the original method (algo 3). This is not memory or computational
% time efficient and it is not reccomended to be used. 
% [L, P, gamma, Nf, eTot] = csspod_original(X, a0, nt0, dt, opts);

%% Plot the eigenspectrum

figure(99)
semilogy(gamma, L(:, 1), 'black-o'); hold on; axis square
semilogy(gamma, L(:, 2), 'blue-o');
semilogy(gamma, L(:, 3), 'red-o');

xlabel('$\gamma$', 'interpreter', 'latex');
ylabel('$\lambda_j$', 'interpreter', 'latex');
legend('$\lambda_1$', '$\lambda_2$', '$\lambda_3$', 'interpreter', 'latex');

%% Plot the CS-SPOD modes (note that these scripts assume that the eigenvectors have been inflated)

% Plot frequency components of the dominant mode
iGam = 10;               % Pick the gamma value to plot
iMode = 1;              % Pick with mode to plot

% Get the current mode to plot
cMode = squeeze(P(iGam, :, :, iMode));

% Frequency components that make up this dominant mode
cF = gamma(iGam) + (-Nf:1:Nf)*a0;

% Plot the real component of the 7-12th frequency components

figure(98)
clf
plot(x, real(cMode(:, 7))); axis square; hold on
plot(x, real(cMode(:, 8)));
plot(x, real(cMode(:, 9))); 
plot(x, real(cMode(:, 10))); 
plot(x, real(cMode(:, 11))); 
plot(x, real(cMode(:, 12))); 
legend({['f = ', num2str(cF(7))],['f = ', num2str(cF(8))], ['f = ', num2str(cF(9))], ['f = ', num2str(cF(10))], ['f = ', num2str(cF(11))], ['f = ', num2str(cF(12))]})
xlabel('x');
ylabel('real(psi(x, f))');
xlim([-50 50]);

% Plot the time dependent modes
disp(['Current gamma value is: ', num2str(gamma(iGam))]);
disp(['Current eigenvalue is: ', num2str(L(iGam, iMode))]);

% Get the current mode to plot
cMode = squeeze(P(iGam, :, :, iMode));
cGam = gamma(iGam);

% Define the time vector you are plotting
tM = 0:dt:(nt0*1-1)*dt; % Plot it for 10 oscillations
pM = tM/(nt0*dt)*2*pi; % Convert to phase

% Get the time domain version of the dominant CS-SPOD mode
pTime = genTimeMode(cMode, cGam, a0, tM, Nf);

% Plot the real and absolute value of the mode for gamma = 0.05
figure(97)
subplot(1, 2, 1);
pcolor(tM, x, real(pTime)); axis square; shading interp; colormap(hot); ylim([-50 50]);
xlabel('t');
ylabel('x');

subplot(1, 2, 2);
pcolor(pM, x, abs(pTime)); axis square; shading interp; colormap(hot); ylim([-50 50]);
xlabel('$\theta$ (rad)');
ylabel('x');
