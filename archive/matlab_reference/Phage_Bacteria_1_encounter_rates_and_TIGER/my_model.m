function dydt = my_model(t,P)

D0 = 2.56;
tauD = 1e-2;
alpha = 1;
K = D0/(tauD^(alpha-1));
D = 6*alpha*K*t^(alpha-1);    %phage diffusivity

vB = 30;
Rb = 1.1;
Ra = 1;


dydt = (5*10^(-5))*[pi*vB*(Rb^2)*P(2) - 4*pi*D*Rb*P(1) - 4*pi*D0*Ra*P(1);
                 -pi*vB*(Rb^2)*P(2) + 4*pi*D*Rb*P(1)];
