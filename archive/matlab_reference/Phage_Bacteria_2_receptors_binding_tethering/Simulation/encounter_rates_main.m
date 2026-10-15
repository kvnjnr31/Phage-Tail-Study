%opengl software
dbstop if error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%simulation parameters%%%%%%%%%%%%%%%%%%%%%%%%%
%sim steps       %start time   %time step
nsteps = 2400;  Time=0;       dt=1e-2;%1e-7;

%simulation space
LX = 15; LY = 15;  %box lengths

x_left=-LX/2; x_right=LX/2; 
y_left=-LX/2; y_right=LX/2; 
z_left=-LX/2; z_right=LX/2; 

%sim volume
V = LX*LX*LX;

%Phage parameters----------------------------------------------------------
%number     capsid radius   tail length     
P = 100;     Rp=0.1;         Ltail = 0.05;                  %lengths in um!                                           
     
%Diffusion parameters
kBT=300*(1.38e-23);                                   %thermal temperature.
eta=.001;                                              %viscocity of water.
Radius = Rp*1e-6;                                             %head radius.
gamma_0 = 6*pi*eta*Radius;                    %coeff. of friction in water.
D0 = (1e12)*kBT./gamma_0;                      %Diffusion constant.

nu = 1;%0.82;                                          %Diffusion Exponent.
tmin = 0.001;                                   %power-law gamma parameter.
m = 0.05*pi; %receptor pitch angle %0.5*bac_length/2*pi; %0.01*pi; %0.14*pi; 
delta = 0.05; %receptor layer thickness
k_bind = 1;
%Bacteria parameters-------------------------------------------------------
%number         %radius     %length             %receptor radius
B = 1000;         Rb=0.8;     bac_length = 1.9;   Rr = 1e-4;
    
%**************************************************************************
%                              Data Metrics
%**************************************************************************
kernel_data = zeros(1,nsteps);                          
mean_head_data = zeros(1,nsteps);
err_head_data = zeros(1,nsteps);
mean_tail_data = zeros(1,nsteps);
err_tail_data = zeros(1,nsteps);
walk_data = zeros(nsteps,3);

%**************************************************************************
%                            Initialize System
%**************************************************************************
%Initialize bacteria receptors
[MatB] = init_receptors(B, delta, m, Rb, bac_length, Rr);

%Initialize phage receptors
[MatP] = init_phage(P, Rb, Rp, x_right, x_left, y_right, y_left, ...
         z_right, z_left, Ltail);
     
%**************************************************************************
%                     Simulate system over nsteps
%**************************************************************************
for n = 1:nsteps

if (mod(nsteps,1000) == 0)
    100*n/nsteps
end

%plot simulation
%simulation_plotter(MatP, MatB, Rb, bac_length, ...
%    kernel_data, mean_head_data, err_head_data, mean_tail_data, err_tail_data, nsteps)

%update phage positions
[MatP] = Langevin_fun(MatP, MatB, Rb, D0, k_bind, bac_length, dt);   
overlap(MatP, Rb, bac_length)

%generate noise and update waiting times    
[MatP] = noise_gen(MatP, D0, dt, tmin, nu, Time, Rb, bac_length);
overlap(MatP, Rb, bac_length)

%Compute Encounters between phage and receptors
[Nenc, Nact, MatP] = binding_fun(MatP, MatB, Rr);

%-------------------------Collect Data Metrics----------------------------
kernel_data(n) = (Nenc/(Nact))*(V/dt);                       

mean_head_data(n) = sum(MatP(:,9),1);
err_head_data(n) = sqrt(var(MatP(:,11),1));

mean_tail_data(n) = median(MatP(:,10),1);
err_tail_data(n) = sqrt(var(MatP(:,10),1));

walk_data(n,:)=MatP(1,1:3);

%apply B.C's---------------------------------------------------------------
[MatP] = boundarys_receptor(MatP, x_right, y_right, z_right, Rp);
overlap(MatP, Rb, bac_length)
%--------------------------------------------------------------------------
Time = Time + dt; %update time

end %time

%
sum(MatP(:,9),1)
%plot simulation
simulation_plotter(MatP, MatB, Rb, bac_length, ...
    kernel_data, mean_head_data, err_head_data, mean_tail_data, err_tail_data, nsteps)

figure(2)
hold on
[X,Y,Z] = cylinder(Rb,1000);
Z(2, :) = bac_length; Z(1, :) = -bac_length;

h = surf(X,Z,Y);
alpha 0.3
set(h,'edgecolor','none','facecolor','b')
hold on

[X,Y,Z] = sphere(1000);
X2 = X(:,1:500) * Rb;  X3 = X(:,500:end) * Rb;
Y2 = Y(:,1:500) * Rb;  Y3 = Y(:,500:end) * Rb;
Z2 = Z(:,1:500) * Rb;  Z3 = Z(:,500:end) * Rb;

h1 = surf(X2,Y2-bac_length,Z2);
alpha 0.3
set(h1,'edgecolor','none','facecolor','b')
h2 = surf(X3,Y3+bac_length,Z3);
alpha 0.3
set(h2,'edgecolor','none','facecolor','b')
hold on

plot3(MatB(:,1),MatB(:,2),MatB(:,3),'ko','MarkerSize',6, 'MarkerFaceColor','k')
plot3(walk_data(10,1),walk_data(10,2),walk_data(10,3),'go','LineWidth',2)
plot3(walk_data(end,1),walk_data(end,2),walk_data(end,3),'ro','LineWidth',2)
plot3(walk_data(10:end,1),walk_data(10:end,2),walk_data(10:end,3),'m','LineWidth',2)
axis(5*[-1 1 -1 1 -1 1])
view(43,24)
%sim_plotter(MatP,MatB,Rb,Rp,LX,LY,Time,N_dep) %plot simulation space
%}