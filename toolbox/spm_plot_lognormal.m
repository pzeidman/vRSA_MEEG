function spm_plot_lognormal(mu_val, v)
% PLOT_NORMAL_LOGNORMAL Plots normal and lognormal distributions
%
% Inputs:
%   mu   - Mean of normal distribution
%   v - Standard deviation of normal distribution

% Define x range for Normal (covering 99% of probability mass)
x  = linspace(-30,30,1000);
L  = spm_Npdf(x,mu_val,v);
V  = spm_Npdf(exp(x),mu_val,v);

fh = figure;
% Normal distribution:
subplot(1,2,1);
area(x,L,'FaceColor',[0.7 0.7 0.7]);
xlabel('lambda');ylabel('Probability density');axis square;
xline(0,':');
set(gca,'FontSize',12);
% Log-Normal distribution:
subplot(1,2,2);
area(exp(x),V,'FaceColor',[0.7 0.7 0.7]);
xlabel('v');ylabel('Probability density');axis square;
xlim([-0.5 10]);
xline(0,':');
set(gca,'FontSize',12);

end
