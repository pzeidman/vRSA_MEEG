function spm_eeg_rsa_review(RSA)

nXt = RSA{1}.M.nXt; % number of basis functions

% Plot parameters
% -------------------------------------------------------------------------

% Ep{bf} = [params x 1]
Ep = cellfun(@(x)spm_vec(x.Ep),RSA,'UniformOutput',false);
Vp = cellfun(@(x)diag(x.Cp),RSA,'UniformOutput',false);

cols = nXt;
rows = 4;

spm_figure('GetWin','vRSA estimated parameters');
spm_clf;

% plot parameters per component
h = [];
i = 1;
for bf = 1:nXt
    h(bf)=subplot(rows,cols,i);
    spm_plot_ci(Ep{bf},diag(Vp{bf}));
    axis square
    title(sprintf('Bf %d',bf));
    set(gca,'FontSize',12);
    if bf == 1
        ylabel('Posterior');
        xlabel('Component');
    end
    i = i + 1;
end
linkaxes(h);

% Plot parameters (exp)
h = [];
for bf = 1:nXt
    h(bf)=subplot(rows,cols,i);
    spm_plot_ci(Ep{bf},diag(Vp{bf}),[],[],'exp');
    axis square
    title(sprintf('Bf %d',bf));
    set(gca,'FontSize',12);
    if bf == 1
        ylabel('Weight');
        xlabel('Component');
    end    
    i = i + 1;
end
linkaxes(h);

% Start new figure
% -------------------------------------------------------------------------
spm_figure('GetWin','vRSA - matrices');
spm_clf;
rows = 3; cols = 4;

% Plot parameter matrix (2nd order)
% -------------------------------------------------------------------------
BB = RSA{1}.BB;
for bf = 2:length(RSA)
    BB = BB + RSA{bf}.BB;
end
subplot(rows,cols,1:2);
imagesc(BB);
axis square
title('Betas (2nd order)','FontSize',16);

% Plot parameters (2nd order)
% -------------------------------------------------------------------------
G = RSA{1}.G;
for bf = 2:length(RSA)
    G = G + RSA{bf}.G;
end
subplot(rows,cols,3:4);
imagesc(G);
axis square
title('vRSA','FontSize',16);

% Plot log evidence for each condition x bf
% -------------------------------------------------------------------------
logBF = cellfun(@(x)spm_vec(x.logBF.cond),RSA,'UniformOutput',false);
logBF = cell2mat(logBF); % params x bfs

nQ = size(logBF,1); % number of components

subplot(rows,cols,6:7);

if nQ == 1
    bar(logBF,'FaceColor','k');
    str = sprintf('Total: %2.2f nats',sum(logBF));
    title({'Log evidence', str},'FontSize',16);
    ylabel('F (with vs without)');
else
    imagesc(logBF);    
    colormap gray
    set(gca,'YTickLabel',RSA{1}.cnames,'YTick',1:nQ);        
    colorbar    
    title({'Log evidence', 'per condition and bf'},'FontSize',16);    
end
xlabel('Basis function');
set(gca,'FontSize',12);
axis square;

% Stop here if there's only one condition-related component
if nQ == 1
    return
end

% Log evidence per condition
% -------------------------------------------------------------------------
F = sum(logBF,2);

subplot(rows,cols,9:10);
bar(F,'FaceColor','k');    
set(gca,'XTickLabel',RSA{1}.cnames(1:nQ),'XTick',1:nQ);    
colormap gray
set(gca,'FontSize',12);
title({'Log evidence', 'per condition'},'FontSize',16);
axis square;

% Log evidence per bf
% -------------------------------------------------------------------------
F = sum(logBF,1);

subplot(rows,cols,11:12);
bar(F,'FaceColor','k');    
xlabel('Basis function');
colormap gray
set(gca,'FontSize',12);
title({'Log evidence', 'per basis function'},'FontSize',16);
axis square;
