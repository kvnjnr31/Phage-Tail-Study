function [] = simulation_plotter(MatP, MatB, Rb, bac_length, ...
    kernel_data, mean_head_data, err_head_data, mean_tail_data, err_tail_data, nsteps)


figure(1)
subplot(2,4,[1 2, 5 6])
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
plot3(MatP(:,1),MatP(:,2),MatP(:,3),'rd','MarkerSize',10, 'MarkerFaceColor','r')
title('Simulation')
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
axis(5*[-1 1 -1 1 -1 1])
%axis([-2.7 2.7 -2.7 2.7 -2.7 2.7])
view(43,24)
hold off
grid off
drawnow

%------------------------------Encounter Kernel----------------------------
step = nsteps/20;
indx = 1:step:nsteps;
t = 1:nsteps;

figure(1)
subplot(2,4,[3 4])
plot(t(indx),kernel_data(indx),'k','LineWidth',2)
title('Encounter Rate')
xlabel('Time')
ylabel('k (cm^3 / hr)')
%axis([0 nsteps 0 1e5])
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow

%----------------------------Head & Tail Statistics------------------------
pp = find(mean_head_data>0);
if ~isempty(pp)

figure(1)
subplot(2,4,[7])
plot(t(indx), (mean_head_data(indx)), 'ko','LineWidth',2) 
title('Head Size')
xlabel('Time')
ylabel('Radius (\mum)')
%axis([0 nsteps -1 1])
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow


figure(1)
subplot(2,4,[8])
plot(t(indx), (mean_tail_data(indx)), 'ko','LineWidth',2)
title('Tail Length')
xlabel('Time')
ylabel('Length (\mum)')
%axis([0 nsteps -1 1])
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow

end






