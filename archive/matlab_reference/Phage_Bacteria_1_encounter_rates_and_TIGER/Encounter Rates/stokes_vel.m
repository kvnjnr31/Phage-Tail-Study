function [S,X,Y] = stokes_vel(v_x,v_y,Rb,LX,LY)
%compute velocity field!!
x = [-LX/2:0.01:LX/2]; y = [-LY/2:0.01:LY/2];
[X, Y] = meshgrid(x,y');

r = sqrt(X.^2+Y.^2);
r(r<Rb)=Rb;
r4 = (r.*r.*r.*r);
dot = v_x*(X) + v_y*(Y);                  %dot product at each pt.(scalar).

v_comp_x = (3*Rb./(4*r)).*(1+(Rb*Rb)./(3*r.*r))*v_x;
v_comp_y = (3*Rb./(4*r)).*(1+(Rb*Rb)./(3*r.*r))*v_y;

x_comp_x = (3*Rb./(4*r)).*(1./(r.*r) - (Rb*Rb)./r4).*(X).*dot;
x_comp_y = (3*Rb./(4*r)).*(1./(r.*r) - (Rb*Rb)./r4).*(Y).*dot;

S = sqrt((v_x-(v_comp_x + x_comp_x)).^2 + (v_y-(v_comp_y + x_comp_y)).^2);

xx = v_x-(v_comp_x + x_comp_x);
yy = v_y-(v_comp_y + x_comp_y);

%contour(X,Y,S,20,'ShowText','on')
[xunitc,yunitc] = circle(0,0,Rb+0.1,50);
[xunit,yunit] = circle(0,0,Rb,1000);
figure(1)
subplot(2,2,[1 3])
colormap('gray');
contourf(X,Y,S,[1 v_x])
hold on
l = streamslice(X,Y,-xx,-yy,0.5);
set(l,'LineWidth',1)
plot(xunit,yunit,'-k','LineWidth',4)
%plot(xunitc,yunitc,'--r','LineWidth',4)
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')

V1 = S((X==0)&(Y>0));
subplot(2,2,2)
plot(y(y>0)-Rb,V1,'k','LineWidth',2)
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold','YAxisLocation','right')
xlabel('Distance (\mum)')
ylabel('v_F (\mum^2/sec)')

V2 = S((X>0)&(Y==0));
subplot(2,2,4)
plot(x(x>0)-Rb,V2,'k','LineWidth',2)
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold','YAxisLocation','right')
xlabel('Distance (\mum)')
ylabel('v_F (\mum^2/sec)')
a==1;
