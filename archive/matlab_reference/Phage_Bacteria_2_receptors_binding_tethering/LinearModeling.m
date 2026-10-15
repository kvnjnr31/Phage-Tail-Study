%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Program to perform linear regression
%                      on adsorption rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Water viscocities (Pa*s)
eta25 = 0.89*(1e-3); 
eta37 = 0.69*(1e-3);

%Phage Diffusion (cm2/hr) T=25C (298K) -> 37C (310K)
DT2 = 2.8*(1e-8)*3600; DT2 = DT2*(310/298)*(eta25/eta37);  %(Zemb)
DLam = 6.2*(1e-8)*3600; DLam = DLam*(310/298)*(eta25/eta37); %(Moldovan)
DT4 = 3.6*(1e-8)*3600;  DT4 = DT4*(310/298)*(eta25/eta37); %(Barr)
DT3 = 1.19*(1e-7)*3600; DT3 = DT3*(310/298)*(eta25/eta37); %(Polson)

D = [DT2 DLam DT4 DT3];

%Adsorption rates (mL/hr)
kT2 = 1.26*(1e-7); %(Puck)
kLam = 9.36*(1e-8); %(Moldovan)
kT4 = 1.38*(1e-7); %(Barr)
kT3 = 1.8*(1e-7); %(Puck)

krates = [kT2 kLam kT4 kT3];

%Bacteria Param's (cm)
Rb = 0.5*(1e-4); %Radius

%Phage hydro radius's (cm)
kbT = 1.38*(1e-23)*298*(3600^2)*(1e4);
eta = eta25*(3600)/(1e2);
Rp = kbT./(6*pi*eta*D);

%Phage tail length's (cm)
rho =[113 150 113 150]*(1e-7);

z = Rb + (Rp +rho);
kmax = 4*pi*D.*(z);
Ps = Rb./(Rb + z);
X = kmax.*Ps;

plot(X,krates,'ko')




