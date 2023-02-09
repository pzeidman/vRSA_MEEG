function [bf,pst] = spm_dcm_eeg2bf(GCM,states)
% Constructs a basis set for convolution by simulating ERPs from a neural
% mass model. The derivative of the ERP is taken with respect to each
% parameter, to form the basis functions.
%
% GCM    - column vector of estimated DCMs
% states - (optional) vector of indices of variables of interest from the
%          DCM
%
% Returns:
% bf  - basis set [time x basis function]
% pst - peri-stimulus time

if nargin < 2
    states = []; % all states
end

Yc = []; % canonical
Yd = []; % derivatives
parfor i = 1:size(GCM,1)
    
    % Unpack DCM
    DCM = GCM{i};    
    M = DCM.M;
    U = DCM.xU;

    % sample parameters from the posterior
    %--------------------------------------------------------------------------
    n_it  = 50;
    Ep = DCM.Ep;
    Cp = DCM.Cp;

    % P: [parameters x iterations]
    rng(1);
    P  = spm_normrnd(spm_vec(Ep),Cp,n_it);
    np = size(P,1);

    % Get time derivative of predicted data (per population)
    %--------------------------------------------------------------------------
    for j = 1:n_it
        fprintf('DCM %d parameterisation %d\n',i,j);
        
        % get parameters
        params = spm_unvec(P, Ep);
        
        % rebuild input structure
        M.fu = @spm_erp_u;
        U.u = feval(M.fu,(1:M.ns)*U.dt,params,M);        

        % get predicted ERP (per latent state and condition)
        % erp{condition} = [time x states]
        erp = spm_gen_erp(Ep,M,U);
        
        % Limit to selected states if requested
        if ~isempty(states)
            for k = 1:length(erp)
                erp{k} = erp{k}(:,states);
            end
        end
        
        % Accumulate canonical responses from all states and conditions
        Yc = [Yc cell2mat(erp')];
        
        % get derivative of ERP (per latent state) w.r.t. parameters
        % dydp{parameter}{condition} = [time x states]
        dydp = spm_diff(@spm_gen_erp,params,M,U,1);
        
        % Limit to selected states if requested
        if ~isempty(states)
            for m = 1:length(dydp)
                for n = 1:length(dydp{m})
                    dydp{m}{n} = dydp{m}{n}(:,states);
                end
            end
        end        

        % Accumulate derivative timeseries from all states and conditions
        for k = 1:np
            Yd = [Yd cell2mat(dydp{k}')];
        end            
        
    end   

end

% Remove empty columns
Yd(:,all(Yd==0)) = [];

% Summarise over parameterisations using SVD (i.e., PCA)
% seperately for canonical and derivative
%--------------------------------------------------------------------------
X = [];

for i = 1:2
    if i == 1
        Y = Yc;
    else
        Y = Yd;
    end
    
    % Mean-correct
    Y = Y - mean(Y);

    % SVD
    [U,S] = spm_svd(Y);
    bf    = U*S;

    % Explained variance by the components
    x = diag(S).^2;
    R = x ./ sum(x);

    % Find the number of components needed to explain 90% of variance
    nc = find(cumsum(R) > 0.9);
    nc = nc(1);
    bf = bf(:,1:nc);

    % Get peristimulus time
    pst = GCM{1}.xY.pst;
    
    X = [X bf];
end

X = spm_orth(X,'norm');
bf = X;