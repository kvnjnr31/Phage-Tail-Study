%load sipho.mat

%B = 1e6/0.01; t = 0.5;
%P = sipho;

%k = -(1/(t*B))*log(1-P(:,2))

load schwartz_data.mat
S = schwartz_data;

rates = S(:,2)*(3600*1e-10);
kmax = max(rates);
Rs = S(:,1); 

%model
Rs0 = abs(rates-0.5*max(rates));
kk = find(Rs0==min(Rs0));
%c = -log((kmax-rates(kk))/kmax)/Rs(kk);
c = 1.25e-2;
%X = [Rs0 ones(length(Rs0),1)];
%Y = -log((kmax - rates0)/kmax);
%beta = inv(X'*X)*(X')*Y;

indx = [1:max(Rs)]';
X = [indx ones(length(indx),1)];
model = kmax*(1-exp(-indx*c));

figure(1)
plot(Rs,rates*1e7,'ko','LineWidth',2)
hold on
plot(indx,model*1e7,'k','LineWidth',2)
xlabel('Receptor Density (receptors / \mum^2)')
ylabel('Adsorbtion Rate (10^{-7}) (mL/hr)')
legend('Data', 'Model')
plot(indx, kmax*1e7+0*indx, 'r--', 'LineWidth',2)
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')




