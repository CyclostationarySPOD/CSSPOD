function P_time = genTimeMode(cMode, cGam, a0, tM, Nf)
% 2nd to end dimensions are spatial
XSize = size(cMode); XSize = XSize(1:end-1);

% Total number of connected frequencies
nf_total = 2*Nf+1;

% Frequency shifts of the data (as a multiple of alpha)
nfM = -Nf:1:Nf;

% Number of times in the time vector
nt_plot = length(tM);

% Flatten input matrix
cMode = reshape(cMode, [prod(XSize), nf_total]);
    
% Current mode in time
P_time = zeros([prod(XSize), nt_plot]);

for p = 1:nf_total
    % Get the current frequency
    cf = cGam + nfM(p)*a0;
    P_time = P_time + exp(1i*2*pi*tM*cf).*cMode(:, p);
end

P_time = reshape(P_time, [XSize, nt_plot]);
end