% ========================================================================
% fig2A_left_repro.m
% Produce a figure verifying bacteria dynamics using:
%   - CTRW simulations for marker values (as in phage_msd.txt)
%   - Empirical MSD curves (dashed) from mucin_diff parameters
% ------------------------------------------------------------------------
% HOW TO RUN:
%   Put this file and phage_msd_ctrw.m in the same folder, then:
%       >> fig2A_left_repro
% ========================================================================

clear; clc; close all; rng(1);

% ---------------- Parameters from mucin_diff -----------------------------
% Mucin concentrations used in calibration (w/v)
Mucin = [1e-6 0.002 0.006 0.01 0.02 0.04];

% Empirical fit parameters for T4 wild-type in mucus
%   MSD_expected(t) = 4 * D_alpha(i) * t^{alpha(i)}
% D_alpha here is in m^2 / s^{alpha} (note the divide-by-4 in the original
% doc to match historical storage).
D_alpha = [15.34 9.5003 4.7175 4.1145 3.5975 1.6845]*(1e-12)/4;
alpha   = [1.0237 0.92924 0.82187 0.81797 0.91178 0.85635];

% The manuscript’s left panel shows ~0%, 0.1%, 1.0%, 4.0%
idx_use = [1 2 4 6];

% Time axis (seconds) used for dashed curves and where we’ll place markers
t_line = logspace(-1, 1, 100);      % smooth line for dashes
t_marks = logspace(-1, 1, 12);      % discrete points for markers

% Conversion to microns^2 for plotting
m2_to_um2 = (1e6)^2;

% ---------------- Figure setup ------------------------------------------
figure('Color','w'); hold on
set(gca,'XScale','log','YScale','log','FontSize',20,...
         'LineWidth',2,'FontWeight','bold')
xlabel('time (sec)')
ylabel('Mean Squared Displacement (\mum^2)')
axis square

% Marker shapes in the same order as the paper
mk = {'o','s','d','^'};   % circle, square, diamond, triangle

% ---------------- Loop over mucin levels --------------------------------
for k = 1:numel(idx_use)
    i = idx_use(k);

    % ----- Expected dashed curve (empirical line) -----
    msd_exp = 4 * D_alpha(i) * (t_line.^alpha(i));      % m^2
    plot(t_line, msd_exp*m2_to_um2, 'k--', 'LineWidth', 2)

    % ----- CTRW simulation for markers (like phage_msd.txt) -----
    %
    % We run a CTRW with waiting-time exponent nu = alpha(i), then
    % amplitude-match the simulated MSD to the expected curve at t_ref
    % (so the markers sit on the calibrated line while keeping CTRW shape).
    %
    nu     = alpha(i);
    dt     = 0.0435;         % as in phage_msd.txt
    steps  = 1e4;            % use 1e4 to get a smooth MSD
    trials = 1e3;            % 1e3 trials (increase if you want smoother)
    [tau_sim, msd_sim] = phage_msd_ctrw(nu, dt, steps, trials);  % m^2

    % Interpolate the simulation to the marker times
    msd_interp = exp(interp1(log(tau_sim), log(msd_sim), log(t_marks), 'linear', 'extrap'));

    % Amplitude-match at a reference time to overlay with empirical curve
    t_ref = 1.0;  % 1 second reference (falls in the center of panel)
    msd_ref_sim = exp(interp1(log(tau_sim), log(msd_sim), log(t_ref), 'linear', 'extrap'));
    msd_ref_exp = 4 * D_alpha(i) * (t_ref^alpha(i));
    scale = msd_ref_exp / msd_ref_sim;
    msd_marks = scale * msd_interp;

    % Plot markers
    plot(t_marks, msd_marks*m2_to_um2, mk{k}, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor','w', ...
        'MarkerSize',8, 'LineWidth',2)
end

% ---------------- Legend (as in the paper) -------------------------------
h = zeros(1,5);
% dummy handles for legend row
for k=1:4
    h(k) = plot(NaN,NaN,mk{k},'MarkerFaceColor','w',...
        'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',8);
end
h(5) = plot(NaN,NaN,'k--','LineWidth',2);
legend(h, {'0%','0.1%','1.0%','4.0%','Expected Value'}, ...
       'Location','northwest','Box','off');

drawnow
