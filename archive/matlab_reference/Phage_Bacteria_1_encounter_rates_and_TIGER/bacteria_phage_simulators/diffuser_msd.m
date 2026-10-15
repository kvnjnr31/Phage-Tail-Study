%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Program to simulate the anomalous diffusion
%                         of phage in mucus via CTRW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parallelzation
matlabpool close force local
clc
NP=12;                                             %number of processors.
matlabpool(NP);

dt = 0.01;                                                    %time step.
sim_length = 5e+3;                          %no. of steps use 1e+3 for MSD.
J = 1e+5;                                %no. of trials. Use 10000 for MSD.

D0 = 1;
alpha = linspace(1e-6,2,20);                          %power-law exponents.

msd_data = zeros(length(alpha),sim_length-1);
exp_data = zeros(length(alpha),1);


parfor kk=1:length(alpha)

alpha_k = alpha(kk); 

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
                tau_c = (dt)*((gamma(1-alpha_k)*U)^(-1/(alpha_k)));
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
%ii = 5:10:(nData-1);
%tau=tau(ii); msd = msd(ii);
x  = log(tau);
y  = log(msd);
p  = polyfit(x,y,1);

nu = p(1);                                              %diffusion exponent.
lnD = p(2);                                    %log - diffusion coefficent.
Dfit = exp(lnD);
yfit = lnD + nu*x;                              

% make logarithmic index
indx0 = round(log(nData));
indx = round(exp(1:0.5:indx0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        plot msd w/ fitted msd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1)
if (kk==2)||(kk==4)||(kk==6)||(kk==8)||(kk==10)
subplot(1,2,1)
plot(log(tau(indx)),log(msd(indx)),'ko','MarkerSize',5,'LineWidth',2)
hold on
plot(log(tau(indx)),log(4*D0/((gamma(alpha_k+1))*(dt^(alpha_k-1)))) ...
    +  alpha_k*log(tau(indx)),'k--','LineWidth',2)
xlabel('log time (sec)')
ylabel('log MSD (\mum^2/s)')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow
end

subplot(1,2,2)
plot(alpha_k,nu,'ko','MarkerSize',5,'LineWidth',2)
hold on
xlabel('Power law exp: \alpha')
ylabel('Diffusion exp: \nu')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow
%}

msd_data(kk,:)=msd;
exp_data(kk,:)=nu;

end

save msd_data.mat msd_data exp_data
matlabpool close
