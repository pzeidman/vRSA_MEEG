function spm_plot_lognormal(mu_val, sigma)
% PLOT_NORMAL_LOGNORMAL Plots normal and lognormal distributions
%
% Inputs:
%   mu_norm   - Mean of normal distribution
%   sigma_norm - Standard deviation of normal distribution
%   mu_log    - Mean of lognormal distribution (log-space)
%   sigma_log - Standard deviation of lognormal distribution (log-space)

% Define x range for Normal (covering 99% of probability mass)
x_norm_min = norminv(0.025, mu_val, sigma);
x_norm_max = norminv(0.975, mu_val, sigma);
x_norm = linspace(x_norm_min, x_norm_max, 1000);

% Compute normal PDF
pdf_norm = normpdf(x_norm, mu_val, sigma);

% Define x range for Lognormal (99% probability mass
x_log_max = logninv(0.9, mu_val, sigma);
x_log = linspace(0, x_log_max, 1000);

% Compute lognormal PDF
pdf_log = lognpdf(x_log, mu_val, sigma);

% Plot Normal Distribution
figure('Name', 'Normal Log normal', 'Units', 'normalized', 'Position', [0.3, 0.3, 0.4, 0.4]); % Adjust figure size
subplot(1,2,1);
fill([x_norm, fliplr(x_norm)], [pdf_norm, zeros(size(pdf_norm))], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(x_norm, pdf_norm, 'k', 'LineWidth', 1);
hold on
plot([mu_val mu_val], [min(pdf_norm), max(pdf_norm)], 'r', 'LineWidth', 2)
xlabel('x');
ylabel('Probability Density');
title(sprintf('N~(mu=%d, sig=%d)', mu_val, sigma));
ylim([min(pdf_norm), max(pdf_norm) + (max(pdf_norm) - min(pdf_norm)) * 0.2])
set(gca, 'box', 'off')
axis square;
set(gca, 'FontSize',12);
% Plot Lognormal Distribution
subplot(1,2,2);
fill([x_log, fliplr(x_log)], [pdf_log, zeros(size(pdf_log))], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(x_log, pdf_log, 'k', 'LineWidth', 1);
ylim([min(pdf_log), max(pdf_log) + (max(pdf_log) - min(pdf_log)) * 0.2])
xlabel('x');
title(sprintf('LN~(mu=%d, sig=%d)', mu_val, sigma));
set(gca, 'box', 'off')
axis square;
set(gca, 'FontSize',12);
set(gcf, 'PaperPositionMode', 'auto'); % Adjust paper size for export
end
