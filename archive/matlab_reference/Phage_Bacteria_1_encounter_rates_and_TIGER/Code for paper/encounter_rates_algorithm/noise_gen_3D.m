%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Program to find jumping phage and generate Brownian noise
%Input parameters
%MatP: Phage matrix
%D0: Diffusion coefficent
%dt: diffusion time
%tmin: minimum waiting time
%nu: Diffusion exponent
%Time: current simulation time

%Output parameters
%%MatP: Phage matrix
%wx, wy, wz: noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MatP, wx, wy, wz] = noise_gen_3D(MatP, D0, dt, tmin, nu, Time)

    kk = find (Time - MatP(:,4)>=0);           %find phage ready to jump.
    
    %make storage for noise
    wx = zeros(size(MatP,1),1); wy = zeros(size(MatP,1),1); wz = zeros(size(MatP,1),1);
    
    %random jump
    Ux = rand(length(kk),1)-0.5; wx(kk) = 0*sqrt(24*D0*dt)*Ux;
    Uy = rand(length(kk),1)-0.5; wy(kk) = 0*sqrt(24*D0*dt)*Uy;
    Uz = rand(length(kk),1)-0.5; wz(kk) = 0*sqrt(24*D0*dt)*Uz;

    %generate waiting times
    U1 = rand(length(kk),1); tau = tmin*((1./U1).^(1./nu));
    MatP(kk,4) = Time + tau;                           %store waiting time.