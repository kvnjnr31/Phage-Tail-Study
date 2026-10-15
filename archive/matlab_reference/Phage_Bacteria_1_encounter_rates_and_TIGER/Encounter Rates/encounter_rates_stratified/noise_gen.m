function [MatP, MatB, wx, wy, wz] = noise_gen(MatP, MatB, dt, ptype, Time, omega,v,thickness, z_right)
%Program to find jumping phage and generate Brownian noise


%find phage ready to jump--------------------------------------------------
    kk = find (Time - MatP(:,4)>=0);    
    mucin = (3/thickness)*(z_right - MatP(kk,3));       %mucin ct at depth.
   [D0,nu,tmin] = diffusion_vs_depth(MatP(kk,3), mucin, ptype);

    
    %make storage for noise
    wx = zeros(size(MatP,1),1); wy = zeros(size(MatP,1),1); wz = zeros(size(MatP,1),1);
    
    %random jump
    Ux = rand(length(kk),1)-0.5; wx(kk) = sqrt(24*D0*dt).*Ux;
    Uy = rand(length(kk),1)-0.5; wy(kk) = sqrt(24*D0*dt).*Uy;
    Uz = rand(length(kk),1)-0.5; wz(kk) = sqrt(24*D0*dt).*Uz;

    %generate waiting times
    U1 = rand(length(kk),1); tau = tmin.*((1./U1).^(1./nu));
    MatP(kk,4) = Time + tau;                           %store waiting time.
    
%update bacteria run direction---------------------------------------------
     
if Time >= MatB(1,7)
    
    theta = 2*pi*rand;                          %random ang: 0 to 2pi.
    alpha = 2*pi*rand;                          %random ang: 0 to 2pi.
    tau_r = -log(rand)/omega;                               %run time.
    MatB(1,7) = Time + tau_r;
    
    MatB(1,4) = abs(v)*cos(alpha)*cos(theta);         %compute vx-velocity.
    MatB(1,5) = abs(v)*cos(alpha).*sin(theta);        %compute vy-velocity.
    MatB(1,6) = abs(v)*sin(alpha);                       %random vz-vector.
    
end 