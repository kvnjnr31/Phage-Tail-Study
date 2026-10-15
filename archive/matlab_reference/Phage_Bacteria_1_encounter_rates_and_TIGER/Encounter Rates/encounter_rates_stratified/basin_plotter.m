function basin_plotter(MatY, MatZ,LX)

%
figure(1)
subplot(2,2,1)
plot(MatY,MatZ,'g.','MarkerSize',5)
ylabel('y')
zlabel('z')
axis([-LX/2 LX/2 -LX/2 LX/2])
axis 'square'
hold on
%pause(0.1)
%}