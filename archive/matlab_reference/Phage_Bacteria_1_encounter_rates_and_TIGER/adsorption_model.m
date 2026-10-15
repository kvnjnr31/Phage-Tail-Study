P10 = 0;
P20 = 10^4;

[t,P] = ode45(@my_model,[0.01 10*60],[P10, P20]); 

figure(1)
hold on
plot(t/60,100*P(:,2)/P20,'k--',t/60,100*P(:,1)/P20,'r--','LineWidth',3)
xlabel('Time: (min)')
ylabel('Free Phage')
axis([0 10 0 105])
set(gca,'FontSize',20,'LineWidth',3,'FontWeight','bold')

figure(2)
hold on
plot(t/60,100*(P(:,2)+P(:,1))/P20,'b--','LineWidth',3)
xlabel('Time: (min)')
ylabel('Free Phage')
axis([0 10 40 105])
set(gca,'FontSize',20,'LineWidth',3,'FontWeight','bold')


