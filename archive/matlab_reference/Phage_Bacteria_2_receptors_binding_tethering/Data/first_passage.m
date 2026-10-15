%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Algorithm to calculate first passage times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tau] = first_passage(kbind,f_rup)

kBT=300*(1.38e-23);                                   %thermal temperature.

%Define model space
x = linspace(0,50,10000)*(1e-9); 
distance = x';
dx = distance(2)-distance(1);

%Define force model
force_model = kbind*distance;                          %define force model.
force_model(force_model<f_rup*(1e-12))=0;            %truncate force model.

%compute potential
dU = -force_model(1:end-1)*dx;
U = cumsum(dU);                                                 %potential.

%Make symmetric landscape
d = [-fliplr(distance(2:end)'),distance(1:end-1)'];              %distance.
d = d';
UU = [fliplr(U'),U'];                                           %potential.
UU = UU';

%Compute escape parameters
Uescape = max(UU)/kBT;                         %energy required for escape.
jj = find(Uescape<=UU/kBT);
xA = min(abs(d(jj)));                         %escape distance (nm).

%--------------Compute First Passage Time Parameters From Model------------
%Diffusion parameters
eta=.001;                                              %viscocity of water.
Radius = 0.1*1e-6;                                            %head radius.
gamma_0 = 6*pi*eta*Radius;                    %coeff. of friction in water.
D = kBT./gamma_0;                                     %Diffusion constant.

%numerically integrate model
U_tilde = UU(abs(d)<=xA)/kBT;            %the potential within xA.

%first integral
dI = exp(U_tilde)*dx;
I = cumsum(dI);
%second integral
dII = exp(-U_tilde).*I*dx;
tau = (1/D)*sum(dII);