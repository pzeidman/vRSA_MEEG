function RSA = spm_eeg_rsa_estimate(RSA)

% Call this function recursively for multiple models
% -------------------------------------------------------------------------
if iscell(RSA)
    if length(RSA) > 1            
        % Recursive call to estimate multiple models
        for i = 1:numel(RSA)
            RSA{i} = spm_eeg_rsa_estimate(RSA{i});            
        end
        
        % Accumulate evidence over models and store in first model
        if size(RSA,2) > 1        
            F = cellfun(@(x)x.logBF.cond',RSA,'UniformOutput',false);
            F = sum(cell2mat(F),2);
            F = F';

            Pp = spm_softmax([F; zeros(1,length(F))]);
            Pp = Pp(1,:);
            
            RSA{1}.allbf.F  = F;
            RSA{1}.allbf.Pp = Pp;
        end
        
        % Get condition related parameters from all bfs [component x bf]
        Ep = cellfun(@(x)full(x.Ep.cond)',RSA,'UniformOutput',false);
        Ep = cell2mat(Ep);
        
        % Generate predicted waveforms and store in first model
        w = sqrt(exp(Ep));        
        RSA{1}.y = RSA{1}.M.Xt * w';

        return
    else
        % Single RSA provided: continue
        RSA = RSA{1};
    end
end
   
% Unpack
% -------------------------------------------------------------------------
X0 = RSA.M.X0;
Xt = RSA.M.Xt;
Q  = RSA.M.Q;
pE = RSA.M.pE;
pC = RSA.M.pC;
Nv = RSA.M.Nv;
BB = RSA.BB;

nq = length(Q);

% Vectorize prior expectations
pE = spm_vec(pE);

% Estimate
% -------------------------------------------------------------------------
[Cy,~,~,F,Fa,Fc,Eh,Ch,hE,hC,Qh] = spm_reml_sc(BB,X0,Q,Nv,pE,pC);

% Scaled matrices as posteriors
Ep = Eh;
Cp = Ch;
Q  = Qh;

% Unscaled matrices as posteriors
%Ep = h;
%Cp = spm_inv(Ph);

% Package results in standard DCM format
% -------------------------------------------------------------------------
% G = 0;
% for i = 1:nq
%     G = G + exp(Ep(i))*Q{i};
% end
G = Cy;

% Unvectorize priors and posteriors
hE = spm_unvec(hE,RSA.M.pE);
Ep = spm_unvec(Ep,RSA.M.pE);

RSA.M.pE = hE;         % prior expectation of parameters
RSA.M.pC = hC;         % prior covariances of parameters
RSA.M.Q  = Q;          % scaled covariance components
RSA.Ep   = Ep;         % posterior expectations
RSA.Cp   = Cp;         % posterior covariance
RSA.G    = G;          % variational similarity matrix
RSA.Fa   = Fa;         % free energy (accuracy only)
RSA.Fc   = Fc;         % free energy (complexity only)
RSA.F    = F;          % free energy

% Supplement with model comparison against null model with no components
% -------------------------------------------------------------------------
rC = RSA.M.pC;
rC(1:end-1,1:end-1)=0;
F0 = spm_log_evidence_reduce(RSA.Ep,RSA.Cp,...
                             RSA.M.pE,RSA.M.pC,...
                             RSA.M.pE,rC);
RSA.F0 = -F0;
    
% Bayesian model comparsion / selection
% -------------------------------------------------------------------------
logBF     = zeros(nq,1);
logBF_BMS = zeros(nq,1);

for q = 1:nq
    
    % Bayesian model comparison (comparison of with/without each component)
    % ---------------------------------------------------------------------
    
    % Switch off just this component by setting a precise shrinkage prior
    rC      = hC;    
    rC(q,q) = 0;
    
    % Get log Bayes factor in favour of reduced model
    % More positive = greater evidence for switching off the parameter
    % More negative = greater evidence for switching on the parameter
    logBF(q) = spm_log_evidence_reduce(RSA.Ep,RSA.Cp,...
                                       RSA.M.pE,RSA.M.pC,...
                                       RSA.M.pE,rC);
    
    % Bayesian model selection (which one component would we retain?)
    % ---------------------------------------------------------------------
    % Switch off all but this component
    rC = hC;
    j = 1:nq;
    j(q) = [];
    rC(j,j) = 0;
    
    % Compare
    logBF_BMS(q,1) = spm_log_evidence_reduce(RSA.Ep,RSA.Cp,...
                                       RSA.M.pE,RSA.M.pC,...
                                       RSA.M.pE,rC);                                       
end

% Flip: more positive = greater evidence for switching on the parameter
logBF = -logBF;

% Log bayes factor -> posterior probability
Pp     = spm_softmax([logBF(:)'; zeros(1,length(logBF))]);
Pp     = Pp(1,:);
Pp_BMS = spm_softmax(logBF_BMS);
    
% Vectors -> structures
RSA.logBF     = spm_unvec(logBF,     RSA.M.pE);
RSA.Pp        = spm_unvec(Pp',       RSA.M.pE);
RSA.Pp_BMS    = spm_unvec(Pp_BMS,    RSA.M.pE);
RSA.logBF_BMS = spm_unvec(logBF_BMS, RSA.M.pE);