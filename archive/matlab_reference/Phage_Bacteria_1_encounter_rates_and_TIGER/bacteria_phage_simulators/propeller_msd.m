%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Algorithm to simulate the movement of bacterial cells.
% Notes: 
% (1) The code will compute the MSD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
%parallelzation
matlabpool close force local
clc
NP=12;                                             %number of processors.
matlabpool(NP);

%system parameters
L1_mat = [1:5:100];                          %array of run rates.
dt = 0.01;%0.01;                                                %time step.
v = 1;                                         %run velocity has mean=24.1.
Sim_Length = 5e+3;                                               %run time.
N = 1e+5;                                                  %no. of trials.

msd_data = zeros(length(L1_mat),Sim_Length-1);
exp_data = zeros(length(L1_mat),1);

parfor k=1:length(L1_mat)
 
    %data storage
    Mat1=zeros(N,2*Sim_Length);                      %storage for location pts.
    100*k/length(L1_mat)
    
    L1 = L1_mat(k);
    
    for n = 1:N
        
    x_data = zeros(1,Sim_Length);                       %storage for plots.
    y_data = zeros(1,Sim_Length);                       %storage for plots.
    
    %initial conditions
    x0b = 0;                                                 %bacteria i.c.
    y0b = 0;
    x_0 = x0b;%((x0b-0.5*v*dt)-(x0b+0.5*v*dt))*rand + (x0b+0.5*v*dt); %random-i.c. 
    y_0 = y0b;%((y0b-0.5*v*dt)-(y0b+0.5*v*dt))*rand + (y0b+0.5*v*dt); %random-i.c.
    
    Time = 0;                                                    %time i.c.
    tau = -log(rand)/L1;                                %transition time i.c.
    T = tau;
    theta = -2*pi*rand + 2*pi;
        
        for j = 1:Sim_Length    

        %bacteria dynamics---------------------------------

           if Time >= T
                
               %select new run direction
                theta = -2*pi*rand + 2*pi;           %random ang: 0 to 2pi.
                tau = -log(rand)/L1;                               %run time.
                T = Time + tau ;
                         
           end
        
        v_x = abs(v)*cos(theta);                       %compute x-velocity.
        v_y = abs(v)*sin(theta);                       %compute y-velocity.

        %advance position    
        x = v_x*dt + x_0;                         
        y = v_y*dt + y_0;            

        x_data(1,j) = x;                         %save x-data for plotting.
        y_data(1,j) = y;                         %save y-data for plotting.
                
        Mat1(n,(2*j-1):(2*j))=[x y];                 %save all experiments.

        x_0 = x;
        y_0 = y;

        Time = Time + dt;

        end
        
    end

% MSD computation of entire ensemble---------------------------------------
nData = Sim_Length;                                    %no. of data points.
msd = zeros(nData-1,1);                                  %storage.
        for h=2:nData 
        particle_new = Mat1(:,1:2);                           %initial pts.  
        particle_old = Mat1(:,(2*h-1):2*h);                   %lag = h*tau.
        deltaCoords = particle_new - particle_old;             %difference.
        squaredDisplacement = sum(deltaCoords.^2,2);            %dx^2+dy^2.
        msd(h-1) = mean(squaredDisplacement(:,1));                 %MSD.
        end
tau0 = 1:(nData-1);
tau = (dt*tau0');

%Linear Fit
iii=floor(0.3*nData):nData-1;


xx = log(tau(iii));
yy = log(msd(iii));
p = polyfit(xx,yy,1);

nu = p(1);
lnD = p(2);

xx2 = dt*(1:(nData-1));
yfit = p(2) + p(1)*log(xx2);                                %p(2)=lnD, p(1)=gamma.

Dfit = exp(lnD)
D_analytical = 2*(v^2)/(L1)
%{
% make logarithmic index
%indx0 = floor(10*log(nData));
%indx = round(exp(0:0.5:indx0/10));
indx=[1:500:nData-1];
figure(1)
subplot(1,2,1)
%if (k == 1)||(k == 2)||(k == 4)||(k == 8)||(k == 20)
if (k == 2)||(k == 3)||(k == 5)||(k == 8)||(k == 15)
plot((tau(indx)), (msd(indx)),'ko','MarkerSize',6)
hold on
plot((tau), D_analytical*tau,'k','LineWidth',2)
end
xlabel('log time (sec)')
ylabel('log MSD (\mum^2/s)')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow

subplot(1,2,2)
if (mod(k,2)~=0)
plot(L1,nu,'ko','MarkerSize',6,'LineWidth',2)
hold on
xlabel('Power law exp: \omega')
ylabel('Diffusion exp: \nu')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow
end
%}

msd_data(k,:)=msd;
exp_data(k,:)=nu;
end

save msd_data.mat msd_data exp_data
matlabpool close
