function plot_fun(kernel, Rcrit, nsteps, v0, D0, dt, nu, omega)

indx = 1;

t=(dt*[1:indx:nsteps-1]);
Dt = 2*((v0^2)/(3*omega))*t + 2*D0*(t.^nu)/(dt^(nu-1));
D = 2*(v0^2)/(3*omega) + 2*nu*D0*(t.^(nu-1))/(dt^(nu-1));
C = 1+(Rcrit^2)./sqrt(pi*Dt.*t);
J = 4*pi*Rcrit*D.*C; %diffusive kernel

Ball = pi*(Rcrit^2)*v0;


figure(1)
hold on
plot(dt*[1:indx:nsteps-1],kernel(1:indx:nsteps-1),'k','LineWidth',3) %mL/hr
%plot(dt*[501:nsteps-1],6*kernel_avg(501:nsteps-1)/(1e+2),'k--','LineWidth',3)  %mL/hr
plot(dt*[1:indx:nsteps-1],J+0*[1:indx:nsteps-1],'b--','LineWidth',3)  %mL/hr
plot(dt*[1:indx:nsteps-1],Ball+0*[1:indx:nsteps-1],'r--','LineWidth',3)  %mL/hr
axis 'square'
xlabel('Time: (sec)')
ylabel('Kernel: \beta (mL/min)')
axis square
%axis([0 10 0 10])
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
