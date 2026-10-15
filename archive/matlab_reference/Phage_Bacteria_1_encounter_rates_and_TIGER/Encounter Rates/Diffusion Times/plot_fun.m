function plot_fun(MatP,Time, D0, Rcrit, nu, dt)

D = (6*D0/(dt^(nu-1)));
t_T4 = ((Rcrit^2)/D)^(1/nu);
%t_T4 = (1/(6*D0))*(Rcrit^2);

n = find(MatP(:,6)==0);
diff_times = mean(MatP(n,5));

figure(1)
hold on
plot([0:dt:Time],t_T4 + 0*[0:dt:Time], 'k','LineWidth',3)
plot(Time,diff_times,'ko','LineWidth',3)
axis 'square'
xlabel('Time: (sec)')
ylabel('Avg. Diffusion Times')
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow
