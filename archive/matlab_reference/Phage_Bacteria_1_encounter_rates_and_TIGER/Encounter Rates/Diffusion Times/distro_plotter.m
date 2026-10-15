opengl software

D0 = 2.54;
t = [0.01:0.01:1000];
R = 1.1;

f = (t./sqrt(4*pi*D0*t)).*exp(-(R^2)./(2*D0*t));
f1 = (1./sqrt(4*pi*D0*t)).*exp(-(R^2)./(2*D0*t));

sum(f)

plot(t,cumsum(f)./cumsum(f1),'k','LineWidth',2)