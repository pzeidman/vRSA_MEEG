function RSA = spm_eeg_rsa_specify(S,D)
% Specifies a variational Representational Similarity Analysis (RSA) model
%
% Overview of vRSA
% -------------------------------------------------------------------------
%
% The EEG/MEG data are modelled according to: Y = X*B + X0*B0 + E, where 
% the rows of the data Y are ordered:
%
%   time 1-T, condition 1
%   time 1-T, condition 2
%   ...
%   time 1-T, condition N
% 
% This function re-expresses the model in terms of the parameters B_hat:
%
%   B_hat = pinv(X)*Y
%
% Where B_hat is the sum of parameters (B), confounds (X0) and noise (E).
% The second order condition-by-condition matrix is then: 
%
%   BB = B_hat * B_hat'
% 
% The objective of vRSA is to decompose matrix BB into experimental effects
% and noise. Contrasts vectors C1,C2,... are specified and then converted to 
% covariance matrices Q1,Q2 according to Qi = Ci * Ci'. The model is then:
%
%   BB = exp(p1)*Q1 + exp(p2)*Q2 + ...
%
% Where p1,p2,... are hyperparameters estimated from the data using
% the function spm_eeg_rsa_estimate.
%
% Inputs
% -------------------------------------------------------------------------
% Required:
% S.Xt          - within-trial design matrix
% S.con_t       - within-trial contrasts (cell array)
% S.con_c       - between-condition contrasts (cell array)
% S.con_t_names - names of within-trial contrasts (optional) 
% S.cnames - names of between-condition contrasts (optional) 
% D             - MEEG data object in standard SPM format
%
% Optional:
% S.X0 - known confounds (optional) 
% S.pE - prior expectation for hyperparameters (optional) 
% S.pC - prior variance for hyperparameters (optional) 
%
% _________________________________________________________________________

% Set defaults
% -------------------------------------------------------------------------

% Required inputs
con = S.con_c;
Xt    = S.Xt;

% Get sizes
[nmodes,ntimes,nconditions] = size(D(:,:,:));
nXt               = size(Xt,2);                % Number of basis functions
RSA.M.nmodes      = nmodes;
RSA.M.ntimes      = ntimes;
RSA.M.nconditions = nconditions;
RSA.M.nXt         = nXt;

% Optional inputs
try cnames = S.cnames; catch, cnames = {}; end
try Xt_names = S.Xt_names; catch, Xt_names = {}; end
try pE          = S.pE;          catch, pE = -16; end
try pV          = S.pV;          catch, pV = 128; end
try X0 = S.X0; catch, X0 = ones(ntimes*nconditions,1); end

if ~iscell(con), con = {con}; end

if isempty(cnames)
    for i = 1:length(con)
        cnames{i} = sprintf('Contrast %d',i);
    end
end

if isempty(Xt_names)
    for i = 1:nXt
        Xt_names{i} = sprintf('BF %d',i);
    end
end

% Prepare the data
% Y: rows encode all times for condition 1, then all times for condition 2
% -------------------------------------------------------------------------

% Y: modes x time x conditions
Y = D(:,:,:);

% Re-shape to Y: time x modes x conditions
Y = permute(Y,[2 1 3]);

% Stack conditions vertically (Y: time*conditions x modes)
splitY = num2cell(Y,[1 2]);
Y = vertcat(splitY{:});

% Mean correct over channels
if size(Y,2) > 1
    Y = Y - mean(Y,2);
else
    Y = Y - mean(Y);
end

% Prepare the design matrix
% -------------------------------------------------------------------------

% X: columns encode design for condition 1, then condition 2, etc
X = kron(eye(nconditions),Xt);

% Identify the basis function of each regressor
bf = repmat((1:nXt)',nconditions,1);

% Estimate first order parameters B
% -------------------------------------------------------------------------
iX = spm_pinv(X);
B  = iX * Y;

% Specify nuisance effects
% -------------------------------------------------------------------------

% Residual projector
XX = [X X0];
R  = speye(size(XX,1)) - XX*spm_pinv(XX);

% Check for rank deficiency
if rank(XX) < size(XX,2)
    warning(['Design matrix is rank deficient. ' ...
             'Ensure regressors don''t sum to a constant']);
end

% estimate spatial degrees of freedom (Nv)
e  = R*Y;
e  = e'*e;
if all(e(:)==0)
    Nv = 1;
else
    Nv = trace(e)^2/trace(e*e);
end

% Projected confounds
if isempty(X0)
    X0 = ones(size(Y,1),1);
end
X0 = iX*X0;

% Prepare contrasts
% -------------------------------------------------------------------------  
% Mean-centre contrasts
con = cellfun(@(C) C - mean(C), con, 'UniformOutput', false);

% Specify covariance components, replicating contrasts over bfs
% -------------------------------------------------------------------------
qnames = {};
qtypes = {};
qcon = [];
qbf = [];

ncon = length(con);
Q = {};
for i = 1:nXt
    for j = 1:ncon
        c = zeros(size(X,2),1);
        c(bf==i) = con{j};
        Q{end+1} = c*c';
        
        qnames{end+1} = sprintf('%s bf(%d)',cnames{j},i);
        qtypes{end+1} = 'Condition';
        qcon(end+1)   = j;
        qbf(end+1)    = i;
    end
end

% Add an extra covariance component for the measurement error
Q{end+1} = iX*iX';
qtypes{end+1} = 'Noise';
qnames{end+1} = 'Noise';
qcon(end+1) = nan;
qbf(end+1)  = nan;

nq = length(Q);

% Set prior covariance matrix
% -------------------------------------------------------------------------
pC = eye(nq) * pV;

% Noise covariance (flat prior)
pC(end,end) = 128;

% Convert prior expectation scalar (pE) to structure
% -------------------------------------------------------------------------
P       = struct();
P.cond  = repmat(pE,1,nq-1);
P.noise = pE;
pE = P;

% Package
RSA.con  = con;        % contrast vectors / matrices
RSA.cnames = cnames;   % names for the contrasts    
RSA.qtypes = qtypes;   % labels for the types of contrast
RSA.qnames = qnames;   % names for the components    
RSA.qcon   = qcon;     % user-specified contrast for each component
RSA.qbf    = qbf;      % basis function for each component
RSA.M.pE = pE;         % prior expectation of parameters
RSA.M.pC = pC;         % prior covariances of parameters
RSA.M.X  = X;          % design matrix
RSA.M.Xt = Xt;         % within-trial design matrix
RSA.M.Xt_names = Xt_names; % Name of the basis functions
RSA.M.bf = bf;         % indices of basis functions within design
RSA.M.X0 = X0;         % null design matrix
RSA.M.Nv = Nv;         % spatial degrees of freedom
RSA.S    = S;          % provided options
RSA.B    = B;          % betas
RSA.Y    = Y;

RSA.BB = RSA.B * RSA.B';

RSA.M.Q = Q;