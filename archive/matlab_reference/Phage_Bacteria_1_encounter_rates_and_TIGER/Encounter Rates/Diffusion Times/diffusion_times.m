opengl software

P = 5e+5;
B = 1;
nsteps = 1000;
Nact = P;

trials = 1;
avg_diff_times=zeros(trials,1);

res_times = [];
%Parameters----------------------------
v=0; Rb=1; Rp=0.1; dt=0.01; Time=0; Rcrit=Rb+Rp;

D0 = 2.54;                           %Diffusion Constant
nu = 1;%0.82;                           %Diffusion Exponent
tmin = dt*sinc(nu)^(1/nu);%0.001;                        %power-law gamma parameter
%***************************************************************************
%                            Main Loop
%***************************************************************************

%***************************************************************************
%                            Initialize System
%***************************************************************************
%initialize phage
for p=1:P 
%store the phage tracks
MatP(p,1) = 0; MatP(p,2) = 0;  MatP(p,3) = 0; %pos-new
%transition timer
MatP(p,4)=0;
%diffusion timer
MatP(p,5)=0;
%active state
MatP(p,6)=1;
end

%***************************************************************************
%                     Simulate system over nsteps
%***************************************************************************
n = 1;

while ((n<nsteps)&&(Nact>0))
    
    %{
    figure(99)
    plot3(MatP(:,1),MatP(:,2),MatP(:,3),'k.')
    axis([-10, 10, -10, 10, -10, 10])
    drawnow
    %}
 
 kk = find (Time - 0*MatP(:,4)>=0);

    for p=1:length(kk) 
       
    %random jump
    Ux = rand-0.5;
    wx = sqrt(24*D0*dt)*Ux;
    Uy = rand-0.5;
    wy = sqrt(24*D0*dt)*Uy;
    Uz = rand-0.5;
    wz = sqrt(24*D0*dt)*Uz;

    MatP(kk(p),1) = MatP(kk(p),1)+wx; MatP(kk(p),2) = MatP(kk(p),2)+wy; MatP(kk(p),3) = MatP(kk(p),3)+wz; 
    end

    %generate waiting times
    U1 = rand(length(kk),1);
    tau = tmin*((1./U1).^(1./nu));
    MatP(kk,4) = Time + tau;                             %store waiting time
 
%Compute Captures--------------------------------------------------------
[MatP] = capture_fun(MatP,Rcrit,Time,dt);
%--------------------------------------------------------------------------
Nact = sum(MatP(:,6))

plot_fun(MatP,Time, D0, Rcrit, nu, dt);

Time = Time + dt; %update time
n = n + 1; %update steps
end %time



