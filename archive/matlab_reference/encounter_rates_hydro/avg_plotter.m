opengl software

%function to plot average residence times
load J_dT4.mat
dT4 = data1;
load J_T4.mat
T4 = data1;

load data_x.mat
x = data2;

figure(1)
subplot(1,2,1)
hold on
plot(x,dT4,'r^','MarkerSize',10,'LineWidth',2)
plot(x,T4,'ko','MarkerSize',10,'LineWidth',2)
axis 'square'
xlabel('Bacteria Velocity: v_b (m/sec)')
ylabel('<Kernel>: \beta (mL/min)')
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow