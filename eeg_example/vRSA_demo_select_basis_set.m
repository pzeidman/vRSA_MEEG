%% Identify the order of Fourier basis set to use

% Get one subject's data
subject = 1;
D = spm_eeg_load(sprintf('subjects/RmD_%s.mat',subjects{subject}));

% Get the first mode of data
t = D.time;
y = squeeze(D(1, :, :));
y = mean(y,2);

nsamples = length(y);

% Orders of Fourier series to try. NB we go one higher than maximum here to
% demonstrate that Matlab's log-likelihood function is working correctly.
orders = 1:ceil(nsamples/2);

% Fit GLMs with fourier basis functions of increasing order
AIC = [];
yhat = {};
for i = 1:length(orders)
    order = orders(i);
    fprintf('Order %d\n',order);
    
    % Create design matrix    
    xBF = struct();
    xBF.dt = 1 / D.fsample;  % Sampling rate
    xBF.name = 'Fourier set';
    xBF.length = nsamples * xBF.dt;
    xBF.order  = order;
    xBF = spm_get_bf(xBF);    
    Xt = xBF.bf(1:nsamples, :);
    
    % Fit GLM. Uses Mathworks Stats Toolbox
    mdl      = fitglm(Xt,y,'linear','Intercept',false);
    yhat{order} = mdl.predict;
    
    % Score for model    
    AIC(order) = mdl.ModelCriterion.AIC;
end

[~,maxorder] = min(AIC);

% Plot
figure;
subplot(1,2,1);
bar(orders,AIC);
xlabel('Fourier basis set order');
title('AIC (lower is better)');
axis square

subplot(1,2,2);
plot(t,y); hold on; plot(t,yhat{maxorder},'--');
xlabel('Time (secs)');
legend({'Data','Prediction'});
title(sprintf('Optimal Fourier order: %d',maxorder));
axis square