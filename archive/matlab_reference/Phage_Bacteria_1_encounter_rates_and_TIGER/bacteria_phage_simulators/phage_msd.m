function [meanD, meanE] = phage_msd(c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Program to simulate the anomalous diffusion
%                         of phage in mucus via CTRW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.0435;                                                    %time step.
sim_length = 1e+4;                          %no. of steps use 1e+3 for MSD.
J = 1e+3;                                %no. of trials. Use 10000 for MSD.

%system parameters
kBT=300*(1.38e-23);                                   %thermal temperature.
eta=.001;                                              %viscocity of water.
Rp = 60e-9;                                                  %phage radius.
gamma_0 = 6*pi*eta*Rp;                        %coeff. of friction in water.
D0 = kBT/gamma_0;

nu = linspace(1e-6,2,20);                             %power-law exponents.

%data storage
%diff_data = zeros(length(cM_k),1);        %storage for diffussion constant.
%exp_data = zeros(length(cM_k),1);         %storage for diffussion exponent.

for kk=1:length(nu)
for k=1:1

nu_k = nu(kk); %1e-50;

Mat1=zeros(J,2*sim_length);                 %storage for particle location.

    for j =1:J
        x_old = 0;                                                  %x-i.c.                   
        y_old = 0;                                                  %y-i.c.
        Time =0;                                                 %time i.c.
        wt=0;                                            %waiting time i.c.
        
        %solve finite diff equation:
        for n=1:sim_length

            if Time >= wt

                wx = rand - 0.5;
                wy = rand - 0.5;
               
                %Langevin x-dir.
                x_new = x_old  + sqrt(24*D0*dt)*wx;
                %Langevin y-dir.
                y_new = y_old  + sqrt(24*D0*dt)*wy;   
                
                x_old=x_new;
                y_old=y_new;
                
                %generate waiting time.  
                U = rand;
                tau_c = (dt)*((gamma(1-nu_k)*U)^(-1/(nu_k)));
                wt = Time + tau_c - dt;                      %waiting time.
           end

        Mat1(j,(2*n-1):(2*n))=[x_old y_old];
        Time = Time + dt;
        
        end
    end

%MSD computation of entire ensemble----------------------------------------
nData = sim_length;                                    %no. of data points.
msd = zeros(nData-1,1);                                  %storage.
        for h=2:nData 
        particle_0 = 0*Mat1(:,1:2);                           %initial pts.  
        particle_h = Mat1(:,(2*h-1):2*h);                   %lag = h*tau.
        deltaCoords = particle_0 - particle_h;             %difference.
        squaredDisplacement = sum(deltaCoords.^2,2);            %dx^2+dy^2.
        msd(h-1) = mean(squaredDisplacement(:,1));    %MSD normalized.
        end
tau = 1:(nData-1);
tau = dt*tau';

%Linear Fit----------------------------------------------------------------
ii = 5:10:(nData-1);
tau=tau(ii); msd = msd(ii);
x  = log(tau);
y  = log(msd);
p  = polyfit(x,y,1);

g = p(1);                                              %diffusion exponent.
lnD = p(2);                                    %log - diffusion coefficent.
Dfit = exp(lnD);
yfit = lnD + g*x;                              

diff_data(k) = exp(lnD)/(2);                      %save diffusion exponent.
exp_data(k) = g;

%--------------------------------------------------------------------------
figure(kk) %plot msd w/ fitted msd
plot(log(tau),log(msd),'ko','MarkerSize',5,'LineWidth',2)
hold on
plot(log(tau),log(4*D0/((gamma(nu_k+1))*(dt^(nu_k-1)))) +  nu_k*log(tau),'k--','LineWidth',2)
%plot(log(tau(ii)),log(4*D0/((gamma(1-nu_k)*gamma(nu_k+1))*(dt^(nu_k-1)))) +  nu_k*log(tau(ii)),'k--','LineWidth',2)
xlabel('lag time (sec)')
ylabel('Ensemble-averaged MSD')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow
hold off
end

D0
Dfit = Dfit/4
g

figure(100)   %plot log-log msd w/ fitted line
plot(nu_k,g,'ko','MarkerSize',10,'LineWidth',2)
hold on
xlabel('Power law exp: \nu')
ylabel('Diffusion exp: \alpha')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow

end
%{
%plot phage track with rings!
ind =1:sim_length;
x_ind = 2*ind - 1;
y_ind = 2*ind;

x_data=Mat1(1,x_ind);
y_data=Mat1(1,y_ind);

figure(3)
plot(x_data,y_data,'k',x_data(1),y_data(1),'go',...
    x_data(end),y_data(end),'ro','LineWidth',1)     %plot run
xlabel('x-position (m)')
ylabel('y-position (m)')
%axis([(x_left) (x_right) (y_lower) (y_upper)])
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow
hold on

for ii=1:3
u1_s=[]; u2_s=[];
[xunit2,yunit2] = circle(0,0,(ii)*(.3e-5),1000);  %draws circles.
u1_s=[u1_s xunit2]; u2_s=[u2_s yunit2];
figure(3)
plot(u1_s,u2_s,'k.','MarkerSize',6)     %plot circles with radii
end
%}