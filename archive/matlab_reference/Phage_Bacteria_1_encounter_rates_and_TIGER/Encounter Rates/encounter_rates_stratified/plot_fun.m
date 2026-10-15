function plot_fun(kernel_new, Rcrit, nsteps, D0, nu, dt)

t=(dt*[20:nsteps-1]).^(nu-1);

%{
%diffusive encounter rates
Dt = D0*t/(dt^(nu-1));
C = 1+(Rcrit^2)./sqrt(pi*Dt.*t);
J = (6/(1e+2))*4*pi*Rcrit*nu*Dt.*C;
%}

%ballistic encounter rates
vb = 30;
J = (6/(1e+2))*pi*(Rcrit^2)*vb;

figure(99)
hold on
plot(dt*[1:20:nsteps-1],6*kernel_new(1:20:nsteps-1)/(1e+2),'k','LineWidth',3) %mL/hr
%plot(dt*[501:nsteps-1],6*kernel_avg(501:nsteps-1)/(1e+2),'k--','LineWidth',3)  %mL/hr
plot(dt*[20:nsteps-1],J+0*[20:nsteps-1],'k--','LineWidth',3)  %mL/hr
axis 'square'
xlabel('Time: (sec)')
ylabel('Kernel: \beta (mL/min)')
axis square
axis([0 10 0 10])
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')


