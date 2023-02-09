function [PEB,F] = spm_eeg_rsa_peb(RSAs,params)
% RSAs - cell array of RSA models [subjects x basis functions]
%np = length(spm_vec(RSAs{1}.Ep));

if nargin < 2
    params = 'all';
end

% Run a PEB model for each RSA
M   = struct();
PEB = spm_dcm_peb(RSAs,M,params);

% Compare PEB models with / without each parameter
F = [];
for i = 1:length(PEB)    
    % conditions x basis functions
    F(:,i) = compare_models(PEB(i));
end

% Accumulate evidence over basis functions and contrasts
F_bfs       = sum(F,1);
F_contrasts = sum(F,2);

Pp_bfs = spm_softmax([F_bfs(:)'; zeros(1,length(F_bfs))]);
Pp_bfs = Pp_bfs(1,:);

Pp_contrasts = spm_softmax([F_contrasts(:)'; zeros(1,length(F_contrasts))]);
Pp_contrasts = Pp_contrasts(1,:);

% Observed and modelled 2nd order
included_bfs = find(Pp_bfs > 0.9);
BB = [];
G  = [];
for i = 1:size(RSAs,1)
    for j = included_bfs
        if i == 1
            BB = RSAs{i,j}.BB;
            G  = RSAs{i,j}.G;
        else
            BB = BB + RSAs{i,j}.BB;
            G  = G + RSAs{i,j}.G;
        end
    end
end

% Stop if requested
if isfield(RSAs{1}.S,'doplot') && RSAs{1}.S.doplot == false
    return
end

% Create plot
spm_figure('GetWin','Group vRSA analysis');
spm_clf;
rows = 3; cols = 2;
ncontrasts = length(RSAs{1}.cnames);

% Plot contrasts
subplot(rows,cols,1);
bar(F_contrasts);
set(gca,'XTick',1:ncontrasts,'XTickLabel',RSAs{1}.cnames);
ylabel('Log evidence (nats)');
set(gca,'FontSize',12);
title('Pooled evidence per component','FontSize',16);

subplot(rows,cols,2);
bar(Pp_contrasts);
set(gca,'XTick',1:ncontrasts,'XTickLabel',RSAs{1}.cnames);
ylabel('Probability');
set(gca,'FontSize',12);
title('Probability per component','FontSize',16);

% Plot basis functions
subplot(rows,cols,3);
bar(F_bfs);
ylabel('Log evidence (nats)');
xlabel('Basis function');
set(gca,'FontSize',12);
title('Pooled evidence per BF','FontSize',16);

subplot(rows,cols,4);
bar(Pp_bfs);
ylabel('Probability');
xlabel('Basis function');
set(gca,'FontSize',12);
title('Probability per BF','FontSize',16);

% Show 2nd order parameters (which serve as data in the RSA)
subplot(rows,cols,5);
%imagesc(BB - diag(diag(BB)));
imagesc(BB);
xlabel('Condition'); ylabel('Condition');
set(gca,'FontSize',12);
title('Betas (second order)','FontSize',16);
axis square;

% Show modelled 2nd order parameters
subplot(rows,cols,6);
%imagesc(G - diag(diag(G)));
imagesc(G);
xlabel('Condition'); ylabel('Condition');
set(gca,'FontSize',12);
title('vRSA','FontSize',16);
axis square;

% -------------------------------------------------------------------------
function [F,P] = compare_models(PEB)
% Family-wise comparisons: each component on vs off
n = length(PEB.Pnames);

% Try switching each parameter off
F = zeros(1,n);
for i = 1:n
    F(i) = reduce_parameters(PEB,i);
end

% Flip log BF to evidence in favour of the full model
F = -F;

% Log BF -> posterior probability for each component individually
P = spm_softmax([F; zeros(1,n)]);
P = P(1,:);

% -------------------------------------------------------------------------
function logBF = reduce_parameters(PEB,q)
% Identifies the evidence for the reduced model with parameter(s) in vector
% q switched off

% Get posteriors
Ep = PEB.Ep;
Cp = PEB.Cp;

% Get priors
pE = PEB.M.pE;
pC = PEB.M.pC;

% Switch off components in q
rC = pC;
rC(q,q) = 0;

% Get log Bayes factor in favour of reduced model
% More positive = greater evidence for switching off the parameter
% More negative = greater evidence for switching on the parameter    
logBF = spm_log_evidence_reduce(Ep,Cp,...
                                pE,pC,...
                                pE,rC);             
