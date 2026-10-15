opengl software
dbstop if error

P = 1e+5; B = 1; nsteps = 1000; trials = 1;

avg_res_times=zeros(trials,1);     %store average resident times per trial. 
avg_kernel=zeros(trials,1);          %store average kernel times per trial.  

data1 = zeros(trials,1);
data2 = zeros(trials,1);

%dimensionality parameters (1=on, 0=off).
D_para_y = 1;                                                       %y-dir.
D_para_z = 1;                                                       %z-dir.  

%sim parameters
omega = 1;                                           %Bacteria Tumble rate.
D0 = 2.54;                                             %Diffusion Constant.
nu = 1;%0.82;                                             %Diffusion Exponent.
tmin = 0.001;                                   %power-law gamma parameter.
dt=0.01; 

v0 = linspace(0,30,trials);                  %array of bacteria velocities.

for indx = 1:trials

indx

res_times = [];              %storage for distributions of residence times. 
dep_vels = [];               %storage for distributions of deposition vels.

%Parameters----------------------------
v=v0(trials); Rb=1; Rp=0.1; td = 0*dt; Time=0; LX = 10; LY = 10;

%simulation space
x_left=-LX/2; x_right=LX/2; 
y_left=-LX/2; y_right=LX/2; 
z_left=-LX/2; z_right=LX/2; 

V = LX*LX*LX; Rcrit=Rb+Rp; x_old =0; y_old =0; z_old =0;                   %Rcrit gives "soft core"

%bacteria runtime
tau = 0;%-log(rand)/omega;                                       %run time.
run_time = Time + tau;

%**************************************************************************
%                            Main Loop
%**************************************************************************
kernel_new = zeros(1,nsteps);                           %initialize kernel.
kernel_avg = zeros(1,nsteps);                   %initialize average kernel.

%**************************************************************************
%                            Initialize System
%**************************************************************************
%--------------------------initialize bacteria-----------------------------
for b=1:B
%fixed bacteria locations
MatB(b,1) = 0; MatB(b,2) = 0; MatB(b,3) = 0; %pos-new
%bacteria velocity
MatB(1,4) = -v; MatB(1,5) = 0; MatB(1,6) = 0;
%bacteria transition timer
MatB(1,7)=0;
end
RR=x_right-0*Rcrit;

%----------------------------initialize phage------------------------------
for p=1:P 
r = 0;

%randomize locations
while r<(Rb+Rp)
x_old=  ((x_right-x_left-2*Rp)*rand + x_left + Rp); %x_right-Rp;
y_old=  ((y_right-y_left-2*Rp)*rand + y_left + Rp); %(2*RR*rand - RR);
z_old=  ((z_right-z_left-2*Rp)*rand + z_left + Rp); %(2*RR*rand - RR);
r = sqrt((x_old-MatB(1,1))^2 + D_para_y*(y_old-MatB(1,2))^2 + D_para_z*(z_old-MatB(1,3))^2);
end

%store the phage tracks
MatP(p,1) = x_old; MatP(p,2) = D_para_y*y_old;  MatP(p,3) = D_para_z*z_old; %pos-new
%transition timer
MatP(p,4)=0;
%activate state
MatP(p,5)=1;
%contact timer
MatP(p,6)=0;

%x-dir old position         %initial x-dir position 
MatP(p,7)=x_old;            MatP(p,10)=x_old;
%y-dir position             %initial y-dir position
MatP(p,8)=y_old;            MatP(p,11)=y_old;
%z-dir position             %initial z-dir position
MatP(p,9)=z_old;            MatP(p,12)=z_old;

end

%Plot histogram of phage concentration distribution------------------------
r = sqrt(MatP(:,1).^2 + MatP(:,2).^2 + MatP(:,3).^2);
%jj = find(MatP(:,1)<0);
%r(jj)=[];

jj = find(r>2.5);
r(jj)=[];
[n,radius] = hist(r,40);

disc_size = (radius(2)-radius(1))/2;
volumes = (4/3)*pi*(radius+disc_size).^3 - (4/3)*pi*(radius-disc_size).^3;

figure(1)
hold on
plot(radius,n./(sum(n)*volumes),'k--','LineWidth',2)
axis square
drawnow
%**************************************************************************
%                     Simulate system over nsteps
%**************************************************************************
n = 1;                                                 %initial step count.
N_dep = 0;

while n<nsteps
%if (mod(n,100)==0)
%   100*n/nsteps 
%end
%{
%jj = find(abs(MatP(:,1))<0.1);
figure(100)
%hold on
%plot3(MatP(jj,1),MatP(jj,2),MatP(jj,3),'k.','MarkerSize',10)
plot(MatP(:,1),MatP(:,2),'k.','MarkerSize',10)
%xlabel('x')
xlabel('y')
ylabel('z')
%axis([-5 5 -5 5 -5 5])
%axis([-4 4 -4 4 -4 4]*0.5)
%view(45, 45)
axis square
pause(0.1)
%}

%update phage positions
[MatP] = Langevin_fun(MatP, MatB, Rb, dt, D_para_y,D_para_z);    
overlap(MatP, MatB, D_para_y, D_para_z, Rb)

%generate noise and update waiting times    
[MatP] = noise_gen(MatP, MatB, D0, dt, tmin, nu, Time, Rb, Rp, D_para_y, D_para_z);
%overlap(MatP, MatB, D_para_y, D_para_z, Rb)

%Compute Encounters
[MatP, MatB, Nenc, Nact, res_times, dep_vels, N_dep] ...
    = encounter_fun_hydro(MatP, MatB, res_times, Rcrit,...
      dep_vels, Rp,Time,D_para_z,D_para_y,dt,LX);
  
kernel_new(n) = (Nenc/(Nact))*(V/dt);                       %update kernel.

if n > 500
kernel_avg(n) = mean(kernel_new(n-500:n));                 %running average.    
end

%update bacteria run direction -> phage position rotation------------------
%select new run direction
     
if 0*Time >= run_time
    
    %generate random unit vector constrained w/i given angle
    angle = 0*3.5; %degrees of change
    vx = (cos(angle*(2*pi/360))-1)*rand + 1;  %x-dir v/random angle
    vz = D_para_z*(2*sqrt((1-vx^2))*rand - sqrt(1-vx^2)); %z-dir
    b = 2*rand - 1;
    vy = (b/abs(b))*sqrt(1 - vx^2 - vz^2); %y-dir
    
    tau = -log(rand)/omega;                                      %run time.
    run_time = Time + tau;
    
    %rotation angles!
    theta = -atan(vy/vx);
    alpha = asin(vz);
    
%rotation about z-axis
x_new     = MatP(:,1)*cos(theta)-MatP(:,2)*sin(theta);
MatP(:,2) = MatP(:,1)*sin(theta)+MatP(:,2)*cos(theta);

%rotation about y-axis
MatP(:,1) =  x_new*cos(alpha)+MatP(:,3)*sin(alpha);
MatP(:,3) = -x_new*sin(alpha)+MatP(:,3)*cos(alpha);

end     
overlap(MatP, MatB, D_para_y, D_para_z, Rb)

%apply B.C's---------------------------------------------------------------
[MatP] = boundarys_hydro(MatP, x_right, y_right, z_right, Rp, D_para_z,D_para_y,RR);
overlap(MatP, MatB, D_para_y, D_para_z, Rb)
%--------------------------------------------------------------------------
Time = Time + dt; %update time
n = n + 1; %update steps

end %time

sim_plotter(MatP,MatB,Rb,Rp,LX,LY,Time,N_dep) %plot simulation space
%plot_fun(kernel_new, Rcrit, nsteps, D0, nu, v, omega, dt)

end

size(MatP)
P

save J_T4.mat data1
save data_x.mat data2

