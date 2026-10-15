%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Program to Perform Linear Modeling
%                             on Moldovan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------Enviromental Parameters----------------------------
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
%Scaling D's
DLam = 6.2*(1e-8)*3600;                    %Phage Diffusion (cm2/hr) T=25C.

%Scaled DLam: T=25C (298K) -> Experimental Temperature (K)
TempK25 = 25+273.15;                                   %Kelvin temp at 25C.
eta25 = 0.89*(1e-3);                               %Water viscosity at 25C.
D = DLam*(T/TempK25).*(eta25./eta);                            %Scaled D's.

%Phage hydro radius's (cm)
kbT = 1.38*(1e-23)*T*(3600^2)*(1e4);
eta = eta*(3600)/(1e2);
Rp = kbT./(6*pi*eta.*D);                                %Phage radius (cm).

rho = 150*(1e-7);                                   %Phage tail length (cm).
a = 0.5*3.8*(1e-4); b = 0.5*0.8*(1e-4); %semimajor and minor axes of bacteria (cm).
Rb = sqrt(a*b/4);                                    %Bacteria Radius (cm).                                   

%Probability of Phage capture
z = Rb + Rp;                                %Starting distance (cm).

%----------------------------Predictors------------------------------------
%Phage Diffusive Encounter Rate
kdiff = 4*pi*D.*z;                 %Diffusive kmax's at experimental temps.

%Bacteria Ballistic Encounter Rate
load('C:\Users\kevin\Documents\MATLAB\VirusProject\Data\MaedaSpeed.mat')
Temp = Maeda(:,1);                                            %temperature.
Vel = Maeda(:,2)*(1e-6);                  %Ecoli velocities at temperature.

Vels = pchip(Temp,Vel,TempsC);                       %splined temperatures.
vb = Vels*(1e2)*3600;                 %(cm/hr) (changes with temperature).
kball = pi*(z.^2).*vb;             %Ballistic kmax's at experimental temps.

%Assign predictor
Pe = (Rb*vb./D);
Sh = 0.5*(1+(1+2*Pe).^(1/3));
%kswim = kdiff.*Sh;                                             %Predictor.
kswim = kdiff;                                                  %Predictor.

%Responses
%----------------Maltose Grown Bacteria (mL / hr)--------------------------
kexp1 = [7.4e-13 6.0e-12 3.6e-12 1.3e-11 3.56e-11 1.4e-11 2.2e-11...
       1.6e-11 2.7e-11 2.6e-11 2.9e-11]*(3600);
   
%----------------Glucose Grown Bacteria (mL / hr)--------------------------
kexp2 = [9.6e-13 3.0e-12 2.3e-12 1.6e-11 2.0e-11 1.2e-11 3.3e-11...
       1.5e-11 2.0e-11 1.6e-11 1.6e-11]*(3600);

%----------------------linear model: kexp1 ~ kmax--------------------------
%parameters
%
outlierx = kswim(5)*(1e7);
outliery = kexp1(5)*(1e7);
kk = kexp1*(1e7); kk(5)=[];
km = kswim*(1e7); km(5)=[];
%}

%{
outlierx = kswim(5);
outliery = kexp1(5);
kk = kexp1; kk(5)=[];
km = kswim; km(5)=[];
%}

tbl = table(kk', km');
tbl.Properties.VariableNames = {'k', 'X'};

%linear fit
mdl = fitlm(tbl,'k ~ X ')
beta = mdl.Coefficients.Estimate;  % the coefficients

figure(1)
x1 = linspace(0, max(km)+1, 1000);
model1 = beta(2)*x1 + beta(1);

[ypred,yci] = predict(mdl,x1','Prediction', 'observation');

subplot(1,2,1)
hold on
plot(x1, model1, 'k--','LineWidth',2,'MarkerSize',10)
plot(x1,yci,'k-.','LineWidth',2)
fill([x1 fliplr(x1)], [yci(:,1)' fliplr(yci(:,2)')], [0.5 0.5 0.5],'facealpha',0.2,'edgecolor','none')
plot(x1,x1,'k', 'LineWidth',2)
plot(km, kk, 'o','MarkerEdgeColor','k', 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10,'LineWidth',2)
plot(outlierx, outliery, 'ko', 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10,'LineWidth',2) 
title('Maltose Grown \it{E. coli}')
xlabel({'Encounter Rate Constant'; '(mL/hr \times 10^{7})'})
ylabel({'Adsorption Rate Constant'; '(mL/hr \times 10^{7})'})
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
%axis([0.4e-7 0.2e-6 0 0.5e-5])
axis([0 2 0 2])
lgd = legend('Data','Outlier','Linear Fit (Maltose)','CI','y = x','Location','northwest');
lgd.FontSize = 12;
axis square
box on

%parameters
%
outlierx = kswim(7)*(1e7);
outliery = kexp2(7)*(1e7);
kk = kexp2*(1e7); kk(7)=[];
km = kswim*(1e7); km(7)=[];
%}

%{
outlierx = kswim(7);
outliery = kexp2(7);
kk = kexp2; kk(7)=[];
km = kswim; km(7)=[];
%}

tbl = table(kk', km');
tbl.Properties.VariableNames = {'k', 'X'};

%----------------------linear model: kexp2 ~ kmax--------------------------
%linear fit
mdl = fitlm(tbl,'k ~ X')
beta = mdl.Coefficients.Estimate;  % the coefficients
%}

figure(1)
x = linspace(0,  max(km)+1, 1000);
model2 = beta(2)*x + beta(1);
[ypred,yci] = predict(mdl,x','Prediction', 'observation');

subplot(1,2,2)
hold on
plot(x, model2, 'k--','LineWidth',2,'MarkerSize',10)
plot(x,yci,'k-.','LineWidth',2)
fill([x fliplr(x)], [yci(:,1)' fliplr(yci(:,2)')], [0.5 0.5 0.5],'facealpha',0.2,'edgecolor','none')
plot(x1,x1,'k', 'LineWidth',2)
plot(km, kk, 'o','MarkerEdgeColor','k', 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10,'LineWidth',2)
plot(outlierx, outliery, 'ko', 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10,'LineWidth',2)
title('Glucose Grown \it{E. coli}')
xlabel({'Encounter Rate Constant'; '(mL/hr \times 10^{7})'})
%ylabel({'Adsorption Rate Constant'; '(mL/hr \times 10^{-7})'})
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
%axis([0.4e-7 0.2e-6 0 0.5e-5])
axis([0 2 0 2])
lgd = legend('Data','Outlier','Linear Fit (Glucose)','CI','y = x','Location', 'northwest');
lgd.FontSize = 12;
axis square
box on

%import phage_bacteria data
Data = readtable('phage_bacteria_data.csv');
A = Data.Adsorption; 
P = Data.Phage;

Adsorption = A(63:71);
Phage = P(63:71);

indx1 = 1; indx2 = 1; indx3 = 1;
%MSM, TSB and NB
for indx = 1:length(Adsorption)
   
    if strcmp(Phage(indx),'T4')
        Phage_T4(indx1) = Adsorption(indx);
        indx1 = indx1 + 1;
    end
    
    if strcmp(Phage(indx),'T2')
        Phage_T2(indx2) = Adsorption(indx);
        indx2 = indx2 + 1;
    end
    
    if strcmp(Phage(indx),'T6')
        Phage_T6(indx3) = Adsorption(indx);
        indx3 = indx3 + 1;
    end
    
end

%T-even phage diffusion constants
D37 = 3.62*3600*(1e-4)^2;                                          %cm2/hr.
TempK24 = 24+273.15;                                   %Kelvin temp at 24C.
TempK37 = 37+273.15;                                   %Kelvin temp at 25C.

D24 = D37*(TempK24/TempK37).*(eta37./eta24);                     %Scaled D.

%T-even phage hydro radius's (cm)
eta = eta24*(3600)/(1e2);
kbT = 1.38*(1e-23)*TempK24*(3600^2)*(1e4);
Rp = kbT./(6*pi*eta*D24);                               %Phage radius (cm).

%Rp = 30*(1e-7);
%D24 = kbT./(6*pi*eta*Rp);

kdiff2 = 4*pi*D24.*(Rb + Rp);

T4Mean = mean(Phage_T4);
T2Mean = mean(Phage_T2);
T6Mean = mean(Phage_T6);

figure(1)
subplot(1,2,1)
plot(kdiff2*ones(1,length(Phage_T4))+(5e-9)*randn(1,length(Phage_T4)), Phage_T4, 'ks', 'MarkerFaceColor','r','MarkerSize',10,'LineWidth',2)
plot(kdiff2*ones(1,length(Phage_T2))+(5e-9)*randn(1,length(Phage_T2)), Phage_T2, 'ks', 'MarkerFaceColor','m','MarkerSize',10,'LineWidth',2)
plot(kdiff2*ones(1,length(Phage_T6))+(5e-9)*randn(1,length(Phage_T6)), Phage_T6, 'ks', 'MarkerFaceColor','b','MarkerSize',10,'LineWidth',2)






