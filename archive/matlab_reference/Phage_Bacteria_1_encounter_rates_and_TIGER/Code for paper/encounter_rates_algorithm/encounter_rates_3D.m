%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Program to Simulate Bac-Phage Encounter Rates
% Uses
%(1) Langevin_fun_3D
%(2) noise_gen_3D
%(3) encounter_fun_3D
%(4) boundarys_3D
%(5) plot_fun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opengl software
dbstop if error

P = 1e+4;                                             %phage concentration.
B = 1;
nsteps = 1e+3;                                            %simulation time.
trials = 1;

%sim parameters

D0 = 2.54;                                             %Diffusion Constant.
nu = 1;%0.82;                                          %Diffusion Exponent.
tmin = 0.001;                                   %power-law gamma parameter.
dt=0.01; 
v0 = 30;                                              %bacteria velocities.
omega = 100; %1                                      %Bacteria Tumble rate.
%Parameters----------------------------
v=v0; Rb=1; Rp=0.1; Time=0; LX = 10; LY = 10;

%simulation space
x_left=-LX/2; x_right=LX/2; 
y_left=-LX/2; y_right=LX/2; 
z_left=-LX/2; z_right=LX/2; 

V = LX*LX*LX; Rcrit=Rb+Rp; x_old =0; y_old =0; z_old =0;

%bacteria runtime
tau = 0;                                                         %run time.
run_time = Time + tau;

%**************************************************************************
%                            Main Loop
%**************************************************************************
kernel = zeros(1,nsteps);                           %initialize kernel.

%**************************************************************************
%                            Initialize System
%**************************************************************************
%initialize bacteria
for b=1:B
%fixed bacteria locations
MatB(b,1) = 0; MatB(b,2) = 0; MatB(b,3) =0; %pos-new
%bacteria velocity
MatB(1,4) = -v; MatB(1,5) = 0; MatB(1,6) = 0;
%bacteria transition timer
MatB(1,7)=0;

%x-dir old position 
MatB(b,8)=MatB(b,1);
%y-dir position 
MatB(b,9)=MatB(b,2);
%z-dir position 
MatB(b,10)=MatB(b,3);
end

%
%initialize phage----------------------------------------------------------
for p=1:P 
r = 0;
while r<(Rb+Rp)
%random locations
x_old= ((x_right-x_left-2*Rp)*rand + x_left + Rp);
y_old= ((y_right-y_left-2*Rp)*rand + y_left + Rp);
z_old= ((z_right-z_left-2*Rp)*rand + z_left + Rp);
r = sqrt((x_old-MatB(1,1))^2 + (y_old-MatB(1,2))^2 + (z_old-MatB(1,3))^2);
end
%store the phage tracks
MatP(p,1) = x_old; MatP(p,2) = y_old;  MatP(p,3) = z_old; %pos-new
%transition timer
MatP(p,4)=0;
%activate state
MatP(p,5)=1;
%contact timer
MatP(p,6)=0;

%x-dir old position 
MatP(p,7)=x_old;
%y-dir position 
MatP(p,8)=y_old;
%z-dir position 
MatP(p,9)=z_old;

end
%--------------------------------------------------------------------------
%}

%**************************************************************************
%                     Simulate system over nsteps
%**************************************************************************
n = 1;
N_dep = 0;

while n<nsteps
    
if (mod(n,100)==0)
   100*n/nsteps 
end

%sim_plotter(MatP,MatB,Rb,Rp,LX,LY,Time,N_dep) %plot simulation space

%generate noise and update waiting times    
[MatP, wx, wy, wz] = noise_gen_3D(MatP, D0, dt, tmin, nu, Time);

%update bacteria - phage positions
[MatP, MatB] = Langevin_fun_3D(MatP, MatB, dt, wx, wy, wz);    

%Compute Encounters
[MatP, MatB, Nenc, Nact] = encounter_fun_3D(MatP, MatB, Rcrit);
kernel(n) = (Nenc/(Nact))*(V/dt);                       %update kernel.

%apply B.C's
[MatP,MatB] = boundarys_3D(MatP, MatB, x_right, y_right, z_right, Rp, Rb);

Time = Time + dt; %update time
n = n + 1; %update steps

end %time

plot_fun(kernel, Rcrit, nsteps, v0, D0, dt, nu, omega)
