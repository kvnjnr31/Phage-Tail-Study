function [tau, msd] = phage_msd_ctrw(nu, dt, sim_length, J)
%PHAGE_MSD_CTRW  CTRW-based anomalous diffusion MSD (ensemble-averaged)
% Mirrors the logic of your phage_msd.txt but returns (tau, msd) without
% plotting, so the caller can place markers where needed.
%
% INPUTS
%   nu          waiting-time power-law exponent (0 < nu <= 2)
%   dt          time step (sec), e.g., 0.0435
%   sim_length  number of steps, e.g., 1e4
%   J           number of trials (ensemble size), e.g., 1e3
%
% OUTPUTS
%   tau         lag times (sec), length sim_length-1
%   msd         ensemble-averaged MSD (m^2) at each tau
%
% Notes:
%   - Uses the same thermal constants seen in the original text file.
%   - We generate isotropic 2D Langevin steps with CTRW waiting times.
%   - MSD is computed against the initial position (as in your code).
%
% Kevin Joiner — NIH prototype reproduction

% physical constants
kBT = 300*(1.38e-23);        % thermal energy
eta = 0.001;                 % viscosity (Pa*s)
Rp  = 60e-9;                 % particle radius (m)
gamma_0 = 6*pi*eta*Rp;       % friction coeff
D0 = kBT/gamma_0;            % baseline diffusion (m^2/s)

% storage
Mat = zeros(J, 2*sim_length);

for j = 1:J
    x = 0; y = 0; Time = 0; wt = 0;
    for n = 1:sim_length
        if Time >= wt
            % Gaussian increments
            wx = rand - 0.5; wy = rand - 0.5;
            x  = x + sqrt(24*D0*dt)*wx;
            y  = y + sqrt(24*D0*dt)*wy;
            % power-law waiting time
            U  = rand;
            tau_c = dt*((gamma(1-nu)*U)^(-1/nu));
            wt = Time + tau_c - dt;
        end
        Mat(j,(2*n-1):(2*n)) = [x y];
        Time = Time + dt;
    end
end

% Ensemble MSD vs lag
nData = sim_length;
msd   = zeros(nData-1,1);
for h = 2:nData
    xy0 = zeros(J,2);
    xyh = Mat(:,(2*h-1):2*h);
    dxy = xyh - xy0;
    msd(h-1) = mean(sum(dxy.^2,2));     % (m^2)
end
tau = (1:(nData-1))' * dt;              % seconds

% (Optional) thinning to reduce points, caller can interpolate as needed
% Keep full vectors here; main script will interpolate to desired markers.
end
