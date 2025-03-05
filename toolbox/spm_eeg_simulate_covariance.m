function Y = spm_eeg_simulate_covariance(Xt, con, s, nmodes, nsub, CV)
% SPM_EEG_SIMULATE_COVARIANCE Simulate time-resolved data with specified covariance components.
%
%   Y = spm_eeg_simulate_covariance(Xt, con, s, nmodes, nsub, CV)
%
%   This function simulates data for n subjects, with covariance components 
%   specified by the contrasts (con) at certain time points, and controlled 
%   by CV. Each subject's response data contains functionally specialized 
%   responses, confounds, and observation error.
%
%   INPUT:
%       Xt      - A time contrast or design matrix (per time point).
%       Encodes the temporal dynamics of the responses to model
%
%       con     - Contrast matrix specifying the relationship between different
%                 conditions or effects to be simulated.
%
%       s       - Standard deviation of the observation noise.
%
%       nmodes  - Number of "modes" or components (e.g., channels, 
%                 principal components...) of the simulated data.
%
%       nsub    - Number of subjects for which data is to be simulated.
%
%       CV      - Covariance indicator, controlling which covariance components
%                 are active at different time points. This can be used to
%                 turn on or off specific contrasts/effects across time.
%
%   OUTPUT:
%       Y       - A cell array of size {nsub, 1}, where each cell contains a
%                 simulated data matrix of size [nTimePoints × nmodes].
%                 The data in Y{i} corresponds to subject i.
%
%   HOW IT WORKS:
%       1. A design matrix (X) is created by taking the Kronecker product 
%          of the contrast (con) with the time-contrast vector/matrix (Xt). 
%          A nuisance effect (an intercept, X0) is also created.
%       2. For each subject:
%          a. A set of functionally specialized responses (B) is generated 
%             from the covariance components (CV) and random noise.
%          b. Observation error (e) is added, scaled by the noise level s.
%          c. Random confounds (B0) are generated for the nuisance intercept X0.
%          d. The final data Y{i} = X*B + X0*B0 + e is constructed and stored.
%
%   EXAMPLE USAGE:
%       % Suppose we have:
%       %   - Xt of size [30 x 2] representing 30 time points and 2 time
%             windows
%       %   - con of size [32 x 2] capturing 2 different contrasts across 32 trials
%       %   - s = 0.1, noise level
%       %   - nmodes = 3
%       %   - nsub = 5 subjects
%       %   - CV specifying which covariance components apply to which time points
%
%       Xt      = randn(30, 2);
%       con     = rand(32, 2) - 0.5;
%       s       = 0.1;
%       nmodes  = 7;
%       nsub    = 10;
%       CV      = [1, 0; 
%                  0, 1];     % Example that toggles the first component in
%                  the first time window and the second component in the
%                  second time window
%
%       Y = spm_eeg_simulate_covariance(Xt, con, s, nmodes, nsub, CV);
%
%       % Y is a 5x1 cell array. Each entry is a 30x3 matrix for one subject.
%
%   See also: SPM, KRON, RANDN
%
%   -----------------------------------------------------------------------

nXt = size(Xt, 2); % N temporal regressors
nC  = size(con, 2); % N contrasts

% Create design matrix and nuisance effect (intercept):
X  = kron(con, Xt); X0 = ones(size(X,1),1); 

% Preallocate Y for each subject
Y  = cell(nsub, 1);

for i = 1:nsub
    % functionally specialised responses, randomly distributed over voxels
    %----------------------------------------------------------------------
    B    = diag(CV(:))*randn(nC * nXt,nmodes);

    % observation error
    %----------------------------------------------------------------------
    e    = randn(size(X,1),nmodes)*s;

    % known confounds
    %----------------------------------------------------------------------
    B0   = randn(size(X0,2),nmodes)/16;

    % response variable
    %----------------------------------------------------------------------
    Y{i} = X*B + X0*B0 + e;
end