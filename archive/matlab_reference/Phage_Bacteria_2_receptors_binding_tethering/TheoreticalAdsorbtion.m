%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Script to Plot Theoretical Adsorption Rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Comparing ballistic v.s. diffusive encounter rates
%--------------------------------------------------------------------------
%Theoretical adsoption rates
Dp = (6.2e-12);                             %Diffusivity of phage (m2/sec).
a = linspace(60*(1e-9),2e-6,1000);
vb = 21*(1e-6);

kdiff = 4*pi*(a)*Dp;                           %Difffusive adsorption rate.
kball = pi*(a.^2)*vb;                            %Ballistic adsorption rate.

kdiff = kdiff*(1e6)*(3600); %(mL / hr)
kball = kball*(1e6)*(3600); %mL / hr)

figure(2)
plot(a,kdiff,'m','LineWidth',2)
hold on
plot(a,kball,'g','LineWidth',2)

%--------------------------------------------------------------------------
%                               Moldovan Paper
%--------------------------------------------------------------------------
%Theoretical adsoption rates
Dp = (6.2e-12);                             %Diffusivity of phage (m2/sec).
a = 0.5*3.8*(1e-6); b = 0.5*0.8*(1e-6); %semimajor and minor axes of bacteria.
b = sqrt(a*b/4);
vb = 30*(1e-6);
L = 150*(1e-9); H = 60*(1e-9);         %phage Lambda tail and head lengths.

kdiff = 4*pi*(b+L+H)*Dp;                       %Difffusive adsorption rate.
kball = pi*((b+L+H)^2)*vb;                      %Ballistic adsorption rate.

kdiff = kdiff*(1e6)*(3600); %(mL / hr)
kball = kball*(1e6)*(3600); %mL / hr)

%Temps of the experiments
Temps = [1 3 4 10 15 20 24 30 35 37 40];

%Maltose Grown Bacteria (mL / hr)
kexp1 = [7.4e-13 6.0e-12 3.6e-12 1.3e-11 3.56e-11 1.4e-11 2.2e-11...
       1.6e-11 2.7e-11 2.6e-11 2.9e-11]*(3600);
   
%Glucose Grown Bacteria (mL / hr)
kexp2 = [9.6e-13 3.0e-12 2.3e-12 1.6e-11 2.0e-11 1.2e-11 3.3e-11...
       1.5e-11 2.0e-11 1.6e-11 1.6e-11]*(3600);

%steric factors
psi1_diff = (kexp1/kdiff)*((b+L+H)/b);
psi2_diff = (kexp2/kdiff)*((b+L+H)/b);

psi1_ball = (kexp1/kball)*((b+L+H)/b);
psi2_ball = (kexp2/kball)*((b+L+H)/b);

z_diff = [psi1_diff psi2_diff];
z_ball = [psi1_ball psi2_ball];

figure(1)
subplot(2,2,[1 3])
plot(Temps, kexp1, 'bo--', 'MarkerEdgeColor','k','MarkerFaceColor','b','LineWidth',2)
hold on
plot(Temps, kexp2, 'ro--', 'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',2)
plot(Temps, kdiff + 0*Temps, 'm','LineWidth',2)
plot(Temps, kball + 0*Temps, 'g','LineWidth',2)
xlabel('Temp (Celsius)')
ylabel('Rate Constant k (mL/hr)')
xt = Temps; xt(2)=[]; xt(10)=[];
xticks(xt)
xticklabels({'1^{\circ}','4^{\circ}','10^{\circ}','15^{\circ}','20^{\circ}','24^{\circ}', '30^{\circ}', '35^{\circ}', '40^{\circ}'})
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')

figure(1)
subplot(2,2,2)
plot(Temps, psi1_diff, 'mo--', 'MarkerEdgeColor','k','MarkerFaceColor','b','LineWidth',2)
hold on
plot(Temps, psi2_diff, 'mo--', 'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',2)
%xlabel('Temp (Celsius)')
ylabel('Steric Factors: \psi')
xt = Temps; xt(2)=[]; xt(10)=[];
%xticks(xt)
%xticklabels({'1^{\circ}','4^{\circ}','10^{\circ}','15^{\circ}','20^{\circ}','24^{\circ}', '30^{\circ}', '35^{\circ}', '40^{\circ}'})
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
axis square

subplot(2,2,4)
plot(Temps, psi1_ball, 'go--', 'MarkerEdgeColor','k','MarkerFaceColor','b','LineWidth',2)
hold on
plot(Temps, psi2_ball, 'go--', 'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',2)
xlabel('Temp (Celsius)')
ylabel('Steric Factors: \psi')
%xt = Temps; xt(2)=[]; xt(10)=[];
%xticks(xt)
%xticklabels({'1^{\circ}','4^{\circ}','10^{\circ}','15^{\circ}','20^{\circ}','24^{\circ}', '30^{\circ}', '35^{\circ}', '40^{\circ}'})
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
axis square

%--------------------------------------------------------------------------
%                               Other Papers
%--------------------------------------------------------------------------

%Experimental Adsorption rates (mL/hr)

kdiff0 = 4*pi*a*Dp/log(2*a/b);                 %Difffusive adsorption rate.
kdiff0 = kdiff0*(1e6)*(3600);

kexp1 = [5.76e-8 9.36e-8 5.94e-7 7.98e-8 1.86e-7 1.26e-7 1.86e-7...
       1.80e-7 1.08e-5 5.40e-7 1.38e-7];
   
%Radius of phage
Rlambda = 30e-9; RT1 = 30e-9; RT2 = 8.77e-7; RT4 = 8.77e-7;
RT3 = 30e-9; RT7 = 30e-9;

%Radii of adsorbed phage
Rp = [Rlambda Rlambda Rlambda Rlambda RT1 RT2 RT4 RT3 RT4 RT4 RT4];

%system parameters
kBT=310*(1.38e-23);                      %thermal temperature K at 37deg C.
%eta=.001;                                             %viscocity of water.
eta = .0006913;                             %viscocity of water at 37deg C.
gamma_0 = 6*pi*eta*Rp;                        %coeff. of friction in water.
D = kBT./gamma_0;                                     %Diffusion of phage.

%Radius of Ecoli
Rb = 1.34e-6;

%Theoretical adsoption rates
kth = (4)*pi*D*Rb*(3600*1e6);

line1 = linspace(0,max(kexp1),100);
plot(kexp1,kth,'ko',line1,line1,'k--')
xlabel('Experimental Adsorption Rate (mL/hr)')
ylabel('Theoretical Adsorption Rate (mL/hr)')


