%Temps of the experiments (C)
TempsC = [1 3 4 10 15 20 24 30 35 37 40];
%Scale Temps to Kelvin (C+273.15 = K)
T = TempsC + 273.15;

%Water viscocities (Pa*s)
eta1 = 1.73*(1e-3); eta3 = 1.62*(1e-3); eta4 = 1.57*(1e-3);
eta10 = 1.31*(1e-3); eta15 = 1.14*(1e-3); eta20 = 1.00*(1e-3);
eta24 = 0.91*(1e-3); eta30 = 0.80*(1e-3); eta35 = 0.72*(1e-3);
eta37 = 0.69*(1e-3); eta40 = 0.65*(1e-3);

eta = [eta1 eta3 eta4 eta10 eta15 eta20 eta24 eta30 eta35 eta37 eta40];

%------------------------Phage-Bacteria param's----------------------------
%Phage hydro radius's (cm)
kbT = 1.38*(1e-23)*T(7)*(3600^2)*(1e4);
eta = eta(7)*(3600)/(1e2);

%Bacteria Parameter
a = 0.5*3.8*(1e-4); b = 0.5*0.8*(1e-4); %semimajor and minor axes of bacteria (cm).
Rb = sqrt(a*b/4);                                    %Bacteria Radius (cm).  

for kk = 1:1
%Phage Radius
Rp =  0.39*(1e-5);                   %Phage Lambda Hydrodynamic Radius (cm).
rho = 150*(1e-7);                   %Phage Lambda typical tail length (cm).

epsilon = linspace(-1e-6,1e-6,1000);              %increase parameter (cm).

Dtail = kbT./(6*pi*eta.*Rp);                                 %Phage Diffusion.
Dhead = kbT./(6*pi*eta.*(Rp+epsilon));                       %Phage Diffusion.                               

%Capture Radius
z = Rb + Rp + rho + epsilon;

%Phage Diffusive Encounter Rate
kdiff_tail = 4*pi*Dtail.*z;                 %Diffusive kmax's at experimental temps.
kdiff_head = 4*pi*Dhead.*z;

figure(1)
hold on
plot(epsilon*(1e7), kdiff_tail*(1e7),'k', epsilon*(1e7), kdiff_head*(1e7),'r','LineWidth',2)
xlabel('\epsilon (nm)')
ylabel('k_{diff} (mL/hr \times 10^{-7})')
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
box on
axis([-12 12 1 2.2])
legend('R_P = constant','\rho = constant','Location', 'northwest')
end

%-------------------------HeatMap!!!---------------------------------------
Dp = linspace(30*(1e-7),200*(1e-7),2000);           %capsid range (cm).
rho = linspace(5*(1e-7), 1000*(1e-7),2000);                %tail range (cm).

[DP,RHO] = meshgrid(Dp,rho);

RP = 0.5*DP;
D = kbT./(6*pi*eta.*RP);                                  %Phage Diffusion.
z = Rb + RP + RHO;
kdiff = 4*pi*D.*z;

figure(2)
hold on
surf(DP*(1e7),RHO*(1e7),(kdiff))
shading interp
xlabel('Capsid Diameter (nm)')
ylabel('Tail Length (nm)')
%zlabel('k_{diff} (mL/hr \times 10^{-7})')
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold','ColorScale','log')
a=colorbar;
a.Limits = [5*(1e-9) 2*(1e-6)];
ylabel(a,'\psi_{diff} (mL/hr)','FontSize',16,'Rotation',270);
a.Label.Position(1) = 5;
view(2)


figure(3)
hold on
c = [.00792    0.0264    0.0120]*(1e-5);
X = DP*(1e7); Y = RHO*(1e7);
[C,h] = contour(X,Y, 0.7*kdiff,c);
h.LineWidth = 2;
h.LineColor = 'k';
xlabel('Capsid Diameter (nm)')
ylabel('Tail Length (nm)')
%zlabel('k_{diff} (mL/hr \times 10^{-7})')
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')

%--------------------------Load Data---------------------------------------
%import Sargasso Sea data
ADSORB = readtable('phage_bacteria_data.csv');
P = ADSORB.Phage;  B = ADSORB.Host; T = ADSORB.Temperature;    

for k = 26 %1:length(P)
    if (strcmp(B{k},'Ecoli')&&(T{k}==24))
        if (strcmp(P{k},'Lambda'))
            Rp =  0.5*63*(1e-7);    %Phage Lambda Hydrodynamic Radius (cm).
            rho = 150*(1e-7);       %Phage Lambda typical tail length (cm).
        end
        
        if (strcmp(P{k},'T1'))
            Rp =  0.39*(1e-5);          %Phage T1 Hydrodynamic Radius (cm).
            rho = 150*(1e-7);           %Phage T1 typical tail length (cm).
        end
        
        if (strcmp(P{k},'T5'))
            Rp =  0.5*65*(1e-5);          %+Phage T5 Hydrodynamic Radius (cm).(Paepe)
            rho = 150*(1e-7);           %Phage T5 typical tail length (cm).
        end
        
        if (strcmp(P{k},'T2'))
            Rp =  27.5*(1e-7);          %Phage T5 Hydrodynamic Radius (cm).
            rho = 290*(1e-7);           %Phage T5 typical tail length (cm).
        end
        
        if (strcmp(P{k},'T4'))
            Rp =  27.5*(1e-7);          %Phage T5 Hydrodynamic Radius (cm).
            rho = 290*(1e-7);           %Phage T5 typical tail length (cm).
        end
        
        if (strcmp(P{k},'T6'))
            Rp =  27.5*(1e-7);          %Phage T5 Hydrodynamic Radius (cm).
            rho = 290*(1e-7);           %Phage T5 typical tail length (cm).
        end
        
        if (strcmp(P{k},'T7'))
            Rp =  0.5*60*(1e-7);        %+Phage T5 Hydrodynamic Radius (cm).(Paepe)
            rho = 290*(1e-7);           %Phage T5 typical tail length (cm).
        end
        
        if (strcmp(P{k},'T3'))
            Rp =  0.5*60*(1e-7);        %+Phage T5 Hydrodynamic Radius (cm).(Paepe)
            rho = 290*(1e-7);           %Phage T5 typical tail length (cm).
        end
        
        if (strcmp(P{k},'J1'))
            Rp =  27.5*(1e-7);          %Phage J1 Hydrodynamic Radius (cm).
            rho = 290*(1e-7);           %Phage J1 typical tail length (cm).
            %Ref: Physicochemical Properties of a Virulent Lactobacillus Phage
            %Containing DNA with Cohesive Ends 
        end
    end
end



