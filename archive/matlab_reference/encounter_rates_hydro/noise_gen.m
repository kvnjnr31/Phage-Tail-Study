function [MatP] = noise_gen(MatP, MatB, D0, dt, tmin, nu, Time, Rb, Rp, D_para_y, D_para_z)

%Program to find jumping phage and generate Brownian noise
    kk = find (Time - MatP(:,4)>=0);             %find phage ready to jump.
    
    %generate waiting times
    U1 = rand(length(kk),1); tau = tmin*((1./U1).^(1./nu));
    MatP(kk,4) = Time + tau;                           %store waiting time.
     
    %random jumps
    for p=1:length(kk) 
    r = 0;

    %randomize locations
    while r<(Rb+Rp)
    Ux = rand-0.5; wx = sqrt(24*D0*dt)*Ux;
    Uy = rand-0.5; wy = sqrt(24*D0*dt)*Uy;
    Uz = rand-0.5; wz = sqrt(24*D0*dt)*Uz;
    
    %update positions
    x_new = MatP(p,1) + wx;                %x-dir
    y_new = D_para_y*(MatP(p,2) + wy);     %y-dir
    z_new = D_para_z*(MatP(p,3) +wz);      %z-dir
    
    r = sqrt((x_new-MatB(1,1))^2 + D_para_y*(y_new-MatB(1,2))^2 + D_para_z*(z_new-MatB(1,3))^2);
    end
    
    MatP(p,1) = x_new; MatP(p,2) = y_new; MatP(p,3) = z_new;
    end

    
    
    