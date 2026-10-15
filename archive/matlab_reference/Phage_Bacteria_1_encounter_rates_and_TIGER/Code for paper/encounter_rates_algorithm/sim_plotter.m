function sim_plotter(MatP,MatB,Rb,Rp,LX,LY,Time,N_dep)

%simplot-------------------------------------------------------------------
%
figure(1)
%subplot(2,2,1)
plot3(MatP(:,1),MatP(:,2),MatP(:,3),'r.','MarkerSize',10)
hold on
plot3(MatB(:,1),MatB(:,2),MatB(:,3),'k.','MarkerSize',20)
xlabel('x')
ylabel('y')
zlabel('z')
view([45,45])
a=LX/2+5;
axis([-a a -a a -a a])
axis 'square'
pause(0.1)
hold off
%}
