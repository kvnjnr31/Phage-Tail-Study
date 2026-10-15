%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Algorithm to simulate the movement of bacterial cells.
% Notes: 
% (1) The code will compute the MSD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
%system parameters
L1_mat = linspace(10,100,20);                          %array of run rates.

dt = 1;                                                      %time step.
v = 24.1e-6;                                   %run velocity has mean=24.1.
Sim_Length = 1e+4;                                               %run time.
N = 1e+4;                                                  %no. of trials.

%data storage
Mat1=zeros(N,2*Sim_Length);                      %storage for location pts.

for ii=1:length(L1_mat)
 
    100*ii/length(L1_mat)
    
    L1 = L1_mat(ii);
    
    for k = 1:N
        
    x_data = zeros(1,Sim_Length);                       %storage for plots.
    y_data = zeros(1,Sim_Length);                       %storage for plots.
    
    %initial conditions
    x0b = 0;                                                 %bacteria i.c.
    y0b = 0;
    x_0 = x0b;%((x0b-0.5*v*dt)-(x0b+0.5*v*dt))*rand + (x0b+0.5*v*dt); %random-i.c. 
    y_0 = y0b;%((y0b-0.5*v*dt)-(y0b+0.5*v*dt))*rand + (y0b+0.5*v*dt); %random-i.c.
    
    Time = 0;                                                    %time i.c.
    T = -dt*log(rand)/L1;                                %transition time i.c.
    theta = -2*pi*rand + 2*pi;
        
        for j = 1:Sim_Length    

        %bacteria dynamics---------------------------------

          % if Time >= T
                
               %select new run direction
                theta = -2*pi*rand + 2*pi;           %random ang: 0 to 2pi.
                n = -log(rand)/L1;                               %run time.
                %T = Time + n ;
                         
        % end
        
        v_x = abs(v)*cos(theta);                       %compute x-velocity.
        v_y = abs(v)*sin(theta);                       %compute y-velocity.

        %advance position    
        x = v_x*n + x_0;                         
        y = v_y*n + y_0;            

        x_data(1,j) = x;                         %save x-data for plotting.
        y_data(1,j) = y;                         %save y-data for plotting.
                
        Mat1(k,(2*j-1):(2*j))=[x y];                 %save all experiments.

        x_0 = x;
        y_0 = y;

        %Time = T;

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
iii=1:nData-1; %floor(0.1*nData)
iiii=1:300:nData-1;

xx = log(tau(iii));
yy = log(msd(iii));
p = polyfit(xx,yy,1);

g = p(1);
lnD = p(2);

yfit = p(2) + p(1)*xx;                                %p(2)=lnD, p(1)=gamma.

Dfit = exp(lnD)
D_analytical = 2*(v^2)/(L1^2)

% make logarithmic index
indx0 = round(log(nData));
indx = round(exp(1:0.5:indx0));

figure(1)

if ii == 1
plot(log(tau(indx)), log(msd(indx)),'ko','MarkerSize',10,'LineWidth',2,'FaceColor','k')
hold on
plot(log(tau(iii)), (yfit),'k--','MarkerSize',2,'LineWidth',2)

elseif ii == 2
plot(log(tau(indx)), log(msd(indx)),'ks','MarkerSize',10,'LineWidth',2)
hold on
plot(log(tau(iii)), (yfit),'k--','MarkerSize',2,'LineWidth',2)

elseif ii == 3
plot(log(tau(indx)), log(msd(indx)),'kd','MarkerSize',10,'LineWidth',2)
hold on
plot((tau(iii)), (yfit),'k--','MarkerSize',2,'LineWidth',2)

elseif ii == 4
plot(log(tau(indx)), log(msd(indx)),'k^','MarkerSize',10,'LineWidth',2)
hold on
plot(log(tau(iii)), (yfit),'k--','MarkerSize',2,'LineWidth',2)

end

xlabel('time (sec)')
ylabel('Mean Squard Displacement \mum^2 sec ^{-1}')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow

end

