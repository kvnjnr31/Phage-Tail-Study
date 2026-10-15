dbstop if error

load binding_force_data.mat
data = binding_force_data;

kBT=300*(1.38e-23);                                   %thermal temperature.
%---------------------Plot Data: Force vs Distance-------------------------
figure(1)
subplot(2,2,[1 3])
hold on
plot([-2:17],-65*[-2:17]-55,'k--','Markersize',12,'Linewidth',2)
plot(data(:,1),data(:,2),'ko','Markersize',8,'Linewidth',1, 'MarkerFaceColor',[0.5 0.5 0.5])
text(20,min(data(:,2)),'1^{st} Rupture ','Color','k','FontSize',14)
text(40,-200,'2^{nd} Rupture ','Color','k','FontSize',14)
text(60,-130,'3^{rd} Rupture ','Color','k','FontSize',14)
title('Binding Force v.s. Distance')
ylabel('Force (pN)')
xlabel('Distance (nm)')
set(gca,'FontSize',16)
legend('Linear Fit','Data','Location','southeast');
axis([-30 90 -1000 50])

%------------------------Compute Full Potential Well-----------------------
distance = (1e-9)*data(:,1); force = (1e-12)*data(:,2);
%Compute potential
dx = (distance(2:end) - distance(1:end-1));
dU = -force(1:end-1).*dx;
U = cumsum(dU);                                                 %potential.
%Make symmetric landscape
d = [-fliplr(distance(1:end-1));distance(1:end-1)];              %distance.
UU = [fliplr(U);U];                                             %potential.

%fit potential to model
sechEqn = 'a*(1-2*exp(b*x)/(exp(2*b*x)+1))';                        %model.
startPoints = [max(UU)/kBT, 0.1];               %initial parameter guesses.
f = fit(d/(1e-9),UU/kBT,sechEqn,'Start', startPoints)

%Compute escape parameters
Uescape = 2635;  %lower 95% confidence interval of model
jj = find(Uescape<UU/kBT);
xA = min(abs(d(jj)/(1e-9)));                              %escape distance.

%Compute diassociation paramters
jj = find(min(force)==force);
Uda = UU(jj)/kBT;                                  %energy at dissociation.
da = d(jj)/(1e-9);                               %distance of dissociation.

%plot model w/data
figure(2)
hold on
legend off
h = plot(f,d/(1e-9),UU/kBT);
%add text to plot
%plot([-xA xA], [Uescape Uescape] ,'b--','Markersize',12,'Linewidth',2)
%plot([-da da], [Uda Uda],'k--','Markersize',12,'Linewidth',2,'MarkerFaceColor','b')
text(20,Uda,'Binding energy at 1^{st} rupture ','Color','k','FontSize',14)
set(h(1), 'MarkerSize',8, 'Marker','o','Color','k','Linewidth',1,'MarkerFaceColor',[0.5 0.5 0.5]); 
set(h(2),'LineWidth',2,'Color','k','LineStyle','--')
title('Receptor Energy Well')
ylabel('Potential Energy: U / k_BT ')
xlabel('Distance (nm)')
set(gca,'FontSize',22)

%extract the handles that require legend entries
hleglines = [h(1) h(2)];
%create the legend
hleg = legend(hleglines,'Data','sech(x) model','Location','southeast');

%------------------------Potential Well Model------------------------------
%Compute potential
kbind = -65*1e-3;
x = linspace(0,50,1000)*1e-9; 
distance = x';
force_model = kbind*distance;
F0 = force_model;
FF0 = [fliplr(F0'),F0'];

dx = distance(2:end) - distance(1:end-1);

%compute the full potential model
dU0 = -force_model(1:end-1).*dx;
U0 = cumsum(dU0);                                   %untruncated potential.
UU0 = [fliplr(U0'),U0'];                                           %potential.
UU0 = UU0';

%truncate force model
f_rup = -79; 
force_model(force_model<f_rup*(1e-12))=0;
FF = [fliplr(force_model'),force_model'];

%compute potential
dU = -force_model(1:end-1).*dx;
U = cumsum(dU);                                                 %potential.

%Make symmetric landscape
d = [-fliplr(distance(2:end)'),distance(1:end-1)'];              %distance.
d = d';
UU = [fliplr(U'),U'];                                           %potential.
UU = UU';

%Compute escape parameters
Uescape = max(UU)/kBT;  %lower 95% confidence interval of model
jj = find(Uescape<=UU/kBT);
xA = min(abs(d(jj)/(1e-9)));                              %escape distance.

%plot model w/data
figure(1)
subplot(2,2,2)
hold on
upp = 2*(1e-9);
low = 0;

indx = (abs(d)<upp)&(abs(d)>low);
plot(d(indx)/(1e-9), UU0(indx)/kBT ,'k--','Markersize',8,'Linewidth',2)
plot(d(indx)/(1e-9), UU(indx)/kBT ,'k','Markersize',8,'Linewidth',2)
%plot(d(indx)/(1e-9), FF0(indx)/(1e-12) ,'k--','Markersize',8,'Linewidth',2)
%plot(d(indx)/(1e-9), FF(indx)/(1e-12) ,'k','Markersize',8,'Linewidth',2)

%add text to plot
txt = ['F_{rupture} = ' num2str(f_rup) ' pN'];
text(2.1,Uescape+1,txt,'Color','k','FontSize',14)
title('Receptor Energy Well Model')
ylabel('U / k_BT ')
xlabel('Distance (nm)')
set(gca,'FontSize',16)
hleg = legend('(1/2)k_{bind} x^2','Truncated','Location','northeastoutside');
title(hleg, 'Models')
%axis([-3 4 -1 32])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Compute Mean First Passage Times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 20;
rupture_forces = linspace(-50,-120,N);   %rupture forces in pico Newtons
tau_data = zeros(1,N);

for kk=1:N
    f_rup = rupture_forces(kk);
    tau_data(kk) = first_passage(kbind,f_rup);
end
 
figure(1)
subplot(2,2,4)
min = 60;
plot(rupture_forces,tau_data/min,'ks','Markersize',8,'Linewidth',2)    
title('Mean Escape Times')
ylabel('<\tau> (min)')
xlabel('F_{rupture} (pN)')
set(gca,'FontSize',16)    

tt = linspace(0,1e-5,1000);
lambda = 1/tau_data(end);
%
figure(3)
hold on
for kk = 1:length(tau_data)
    
    lambda = 1/tau_data(kk);
    plot(tt,lambda*exp(-lambda.*tt),'k','Markersize',8,'Linewidth',2)
    
end
title('Mean Escape Times')
ylabel('Mean Escape Time: <\tau> (sec)')
xlabel('F_{rupture}')
set(gca,'FontSize',22)      
    
    


