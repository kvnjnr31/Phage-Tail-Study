%------------------Phage Paramters v.s temperature-------------------------
%Temps of the experiments (C)
TempsC1 = [1 3 4 10 15 20 24 30 35 37 40];
%Scale Temps to Kelvin (C+273.15 = K)
T = TempsC1 + 273.15;

%Water viscocities (Pa*s)
eta1 = 1.73*(1e-3); eta3 = 1.62*(1e-3); eta4 = 1.57*(1e-3);
eta10 = 1.31*(1e-3); eta15 = 1.14*(1e-3); eta20 = 1.00*(1e-3);
eta24 = 0.91*(1e-3); eta30 = 0.80*(1e-3); eta35 = 0.72*(1e-3);
eta37 = 0.69*(1e-3); eta40 = 0.65*(1e-3);

eta = [eta1 eta3 eta4 eta10 eta15 eta20 eta24 eta30 eta35 eta37 eta40];

%Scaling D's
DLam = 6.2*(1e-8)*3600;                    %Phage Diffusion (cm2/hr) T=25C.

%Scaled DLam: T=25C (298K) -> Experimental Temperature (K)
TempK25 = 25+273.15;                                   %Kelvin temp at 25C.
eta25 = 0.89*(1e-3);                               %Water viscosity at 25C.
DLams = DLam*(T/TempK25).*(eta25./eta);                        %Scaled D's.

%Compute phage hydro radius's (cm)
kbT = 1.38*(1e-23)*T*(3600^2)*(1e4);
eta = eta*(3600)/(1e2);
Rp = kbT./(6*pi*eta.*DLams);
Rp = Rp(1);

%Phage-Bacteria param's
rho =150*(1e-7);                                   %Phage tail length (cm).
a = 0.5*3.8*(1e-4); b = 0.5*0.8*(1e-4); %semimajor and minor axes of bacteria (cm).
Rb = sqrt(a*b/4);                                    %Bacteria Radius (cm).
                                   
%Probability of Phage capture
z = Rb + (Rp +rho);                                %Starting distance (cm).

%------------------Bacteria Parameters vs temperature----------------------
load('C:\Users\kevin\Documents\MATLAB\VirusProject\Data\MaedaSpeed.mat')
Temp = Maeda(:,1);
Vel = Maeda(:,2)*(1e-6);

load('C:\Users\kevin\Documents\MATLAB\VirusProject\Data\MaedaTumble.mat')
Tumb = MaedaTumble(:,2);

%Temps of the experiments (C)
TempsC = [1 3 4 10 15 20 24 30 35 37 39]; %linspace(1,39.1376,100); %[1 3 4 10 15 20 24 30 35 37 40];
Vels = pchip(Temp,Vel,TempsC);                         %splined velocities.
Tumbs = pchip(Temp,Tumb,TempsC);                          %splined tumbles.

%Bacteria Diffusion(s)
vb = Vel*(1e2); %(cm/sec) (changes with temperature)
tau = 1./Tumb;

vbs = Vels*(1e2); %(cm/sec) (changes with temperature)
taus = 1./Tumbs;

Decoli = (1/3)*(tau.*vb.^2)*3600;      %(cm2/hr)
Decolis = (1/3)*(taus.*vbs.^2)*3600;   %(cm2/hr)

figure(1)
subplot(1,2,1)
plot(Temp, Vel*(1e6),'o','MarkerEdgeColor','k', 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10,'LineWidth',2)
hold on
plot(TempsC, Vels*(1e6), 'ko--','MarkerFaceColor','k','LineWidth',2,'MarkerSize',10)
xlabel('Temperature deg C')
ylabel('v_{E. coli} (\mum/sec)')
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
legend('Data','Spline','Location', 'northwest')
axis square

subplot(1,2,2)
semilogy(TempsC1, DLams, 'o','MarkerEdgeColor','k', 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10,'LineWidth',2)
xlabel('Temperature deg C')
ylabel('D_{phage \lambda}(cm^2/hr)')
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
axis([0 40 0.9e-5 1e-3])
axis square

figure(2)
subplot(1,2,1)
semilogy(Temp, Decoli, 'ko','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2)
hold on
semilogy(TempsC, Decolis, 'k--','MarkerEdgeColor','k','LineWidth',2)
semilogy(TempsC1, DLams, 'o--','MarkerEdgeColor','k', 'MarkerFaceColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5],'LineWidth',2)
xlabel('Temperature deg C')
ylabel('D (cm^2/hr)')
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
axis([0 40 1e-5 1])
axis square

figure(2)
subplot(1,2,2)
semilogy(TempsC, 4*3*pi*(z).*Decolis, 'b--','MarkerEdgeColor','k','LineWidth',2)
hold on
semilogy(Temp, 4*3*pi*(z).*Decoli, 'ko','MarkerEdgeColor','k','MarkerFaceColor','b','LineWidth',2)
semilogy(TempsC1, 4*pi*(z).*DLams, 'ko--','MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5],'LineWidth',2)
semilogy(TempsC, pi*((z).^2).*vbs*3600, 'r--','LineWidth',2)
semilogy(Temp, pi*((z).^2).*vb*3600, 'ro','MarkerFaceColor','r','LineWidth',2)
xlabel('Temperature deg C')
ylabel('k_{max} (cm^3/hr)')
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
axis([0 40 1e-10 1e-1])
axis square


%--------------Typical length of bacteria trajectories---------------------
d = Vel.*tau;
figure(2)
plot(Temp,d./(1e-7),'ko-')
xlabel('Temperature deg C')
ylabel('Phage per E.Coli run length')
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')


%Temps of the experiments (C)
TempsC = [1 3 4 10 15 20 24 30 35 37 40];

%Water viscocities (Pa*s) (1cp = mPa*s)
eta1 = 1.73*(1e-3); eta3 = 1.62*(1e-3); eta4 = 1.57*(1e-3);
eta10 = 1.31*(1e-3); eta15 = 1.14*(1e-3); eta20 = 1.00*(1e-3);
eta24 = 0.91*(1e-3); eta30 = 0.80*(1e-3); eta35 = 0.72*(1e-3);
eta37 = 0.69*(1e-3); eta40 = 0.65*(1e-3);

eta = [eta1 eta3 eta4 eta10 eta15 eta20 eta24 eta30 eta35 eta37 eta40];

figure(3)
plot(TempsC, eta*(1e3), 'ko-','MarkerEdgeColor','k','LineWidth',2)
xlabel('Temperature deg C')
ylabel('Viscosity (cP)')
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')

%Phage Diffusion (cm2/hr) T=25C. (changes with viscocity)
DLam = 6.2*(1e-8)*3600;                                      %Phage Lambda.


