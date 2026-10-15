function plot_fun(kernel_new, Rcrit, nsteps, D0, nu, v, omega, dt)

%Plot Encounter Rates----------------------------------------------
t=(dt*[1:20:nsteps-1]);
Dt = 2*((v^2)/(3*omega))*t + 2*D0*(t.^nu)/(dt^(nu-1));
D = 2*(v^2)/(3*omega) + 2*nu*D0*(t.^(nu-1))/(dt^(nu-1));
C = 1+(Rcrit^2)./sqrt(pi*Dt.*t);

J = 4*pi*Rcrit*D.*C;                                     %diffusive kernel.
Ball = pi*(Rcrit^2)*v;                                   %ballistic kernel.

subplot(1,2,2)
hold on
plot(dt*[1:20:nsteps-1],6*kernel_new(1:20:nsteps-1)/(1e+2),'b','LineWidth',3) %mL/hr
%plot(dt*[1:20:nsteps-1],6*J/(1e+2)+0*[1:20:nsteps-1],'b','LineWidth',3)  %mL/hr
plot(dt*[1:20:nsteps-1],6*Ball/(1e+2)+0*[1:20:nsteps-1],'r','LineWidth',3)  %mL/hr
axis 'square'
xlabel('Time: (sec)')
ylabel('Kernel: \beta (mL/min)')
axis square
axis([0 10 0 50])
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')


%{
%Generate PDF----------------------------------------------
d=find(dep_vels<0);
[nout1 nbins1]=hist(res_times,max(res_times)/dt);
[nout2 nbins2]=hist((dep_vels),50);
figure(2)
subplot(1,2,1)
hold on
plot(nbins1, nout1/((nbins1(2)-nbins1(1))*sum(nout1)), 'ko','LineWidth',2)
xlabel('T_{res}')
ylabel('PDF')
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow

subplot(1,2,2)
hold on
plot(nbins2, nout2/((nbins2(2)-nbins2(1))*sum(nout2)), 'ko','LineWidth',2)
xlabel('v_p^n')
ylabel('PDF')
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow

figure(3)
hold on
plot((res_times), (dep_vels), 'ko','LineWidth',3)
xlabel('T_{res}')
ylabel('-v_p^n')
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
axis square
drawnow
%}





