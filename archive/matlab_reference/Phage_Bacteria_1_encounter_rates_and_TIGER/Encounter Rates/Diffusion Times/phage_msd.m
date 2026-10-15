%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Program to simulate the anomalous diffusion
%                         of phage in mucus via CTRW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.01;                                                    %time step.
sim_length = 1e+4;                          %no. of steps use 1e+3 for MSD.
J = 1e+3;                                %no. of trials. Use 10000 for MSD.

%system parameters
D0 = 2.54;
nu = 0.82;                                             %power-law exponent.
tmin = dt*sinc(nu)^(1/nu);
Mat1=zeros(J,3*sim_length);                 %storage for particle location.

    for j =1:J
        x_old = 0;                                                  %x-i.c.                   
        y_old = 0;                                                  %y-i.c.
        z_old = 0;                                                  %z-i.c.
        Time =0;                                                 %time i.c.
        wt=0;                                            %waiting time i.c.
        
        %solve finite diff equation:
        for n=1:sim_length

            if Time >= wt

                wx = rand - 0.5;
                wy = rand - 0.5;
                wz = rand - 0.5;
               
                %Langevin x-dir.
                x_new = x_old  + sqrt(24*D0*dt)*wx;
                %Langevin y-dir.
                y_new = y_old  + sqrt(24*D0*dt)*wy; 
                %Langevin z-dir.
                z_new = z_old  + sqrt(24*D0*dt)*wz;
                
                x_old=x_new;
                y_old=y_new;
                z_old=z_new;
                
                %generate waiting time.  
                U = rand;
                tau_c = (tmin)*(1/U)^(1/(nu));
                wt = Time + tau_c;                          %waiting time.
           end

        Mat1(j,(3*n-2):(3*n))=[x_old y_old z_old];
        Time = Time + dt;
        
        end
    end

%MSD computation of entire ensemble----------------------------------------
nData = sim_length;                                    %no. of data points.
msd = zeros(nData-1,1);                                           %storage.
        for h=2:nData 
        particle_0 = 0*Mat1(:,1:3);                           %initial pts.  
        particle_h = Mat1(:,(3*h-2):3*h);                   %lag = h*tau.
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

%diff_data(k) = exp(lnD)/(2);                      %save diffusion exponent.
%exp_data(k) = g;

%--------------------------------------------------------------------------
figure(1) %plot msd w/ fitted msd
plot(log(tau),log(msd),'ko','MarkerSize',5,'LineWidth',2)
hold on
plot(log(tau),log(6*D0/(dt^(nu-1))) +  nu*log(tau),'k--','LineWidth',2)
xlabel('lag time (sec)')
ylabel('Ensemble-averaged MSD')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow
hold off

figure(2) %plot msd w/ fitted msd
plot((tau),(msd),'ko','MarkerSize',5,'LineWidth',2)
hold on
plot((tau),(6*D0/(dt^(nu-1))).*(tau.^nu),'k--','LineWidth',2)
xlabel('lag time (sec)')
ylabel('Ensemble-averaged MSD')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow
hold off

%{
figure(100)   %plot log-log msd w/ fitted line
plot(nu_k,g,'ko','MarkerSize',10,'LineWidth',2)
hold on
xlabel('Power law exp: \nu')
ylabel('Diffusion exp: \alpha')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow
%}

