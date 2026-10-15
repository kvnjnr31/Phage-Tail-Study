%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Program to compute Sherwood number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opengl software
dbstop if error

Rc = 1.1; vb = 30; tau = 1/100;
Db = (tau*vb^2)/3;
Dt4 = 0; Ddt4 = 2.54; 
D_phage = Ddt4;
kt4 = 2.8e-7; kdt4 = 1.4e-7; bt4 = 1.49e-7; bdt4 = 2.13e-7;

J=4*pi*(0*Db + D_phage)*Rc*(60*60*(1e6))*1e-18

U=linspace(0,50,100);  %parcle velocity
L = 1;                 %particle diameter

Dwt = [3.84, 2.38, 1.18, 1.03, 0.90, 0.42]; %diffusion coeff's
Dhoc = [3.64, 3.14, 2.63, 2.54, 1.37, 0.74]; %diffusion coeff's

for k=4
Pe_wt = U*L./Dwt(k);
Sh_wt = 0.5*(1+(1+Pe_wt).^(1/3)); %Sherwood number
J_wt = 4*pi*Dwt(k)*Sh_wt;

Pe_hoc = U*L./Dhoc(k);
Sh_hoc = 0.5*(1+(1+Pe_hoc).^(1/3)); %Sherwood number
J_hoc = 4*pi*Dhoc(k)*Sh_hoc;

gray = [0.5 0.5 0.5];
figure(1)
subplot(1,2,1)
hold on
plot(U,Sh_wt,'k','LineWidth',2)
plot(U,Sh_hoc,'Color',gray,'LineWidth',2)
xlabel('Relative Swimming Speed v/d')
ylabel('Sherwood Number')
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
%axis([0 30 0 3])
axis square;

subplot(1,2,2)
hold on
plot(U,J_wt,'k','LineWidth',2)
plot(U,J_hoc,'Color',gray,'LineWidth',2)
xlabel('Relative Swimming Speed v/d')
ylabel('Encounter Kernel: \beta')
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
axis([0 30 0 100])
axis square;
end