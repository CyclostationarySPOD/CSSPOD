function opts = defaultparams(X, nt0)
% This function generates the default parameters. We advice the user to use
% their own parameters.

disp('Setting default parameters required for CS-SPOD');

% Size of X
XS                  = size(X);

% Hamming window with 20*nt0 is the default window
opts.Window         = hamming(nt0*20);

% Default: unit weight matrix 
% Must be a matrix of (nx1, nx2, nx3, 1) size
opts.Weight         = ones([XS(2:end), 1]);

% Default: block overlap
opts.FracOvlp       = 0.67;

% Default: use double sided spectrum
opts.singlef        = 0;

% Default: save 5 modes
opts.nsave          = 5;

% Default: do not subtract block mean
opts.blockwisemean  = false;

% Default: use as many connected frequencies as the data allows
opts.Nf             = 'max';

% Default: display output messages
opts.quiet          = 0;

% Default: inflate the output eigenvectors
opts.inflate        = 1;

end
