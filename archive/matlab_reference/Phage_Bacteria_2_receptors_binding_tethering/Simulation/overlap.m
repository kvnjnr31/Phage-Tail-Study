%check for overlaps

function [] = overlap(MatP, Rb, bac_length)

for p=1:size(MatP,1) 
r_body = sqrt((MatP(p,1))^2 + (MatP(p,3))^2);
r_head = sqrt((MatP(p,1))^2 + (MatP(p,2)-bac_length)^2 + (MatP(p,3))^2);
r_tail = sqrt((MatP(p,1))^2 + (MatP(p,2)+bac_length)^2 + (MatP(p,3))^2);

%randomize locations
if ((r_body<Rb)&&(abs(MatP(p,2))<bac_length))||((r_head<Rb)&&(MatP(p,2)>bac_length))||((r_tail<Rb)&&(MatP(p,2)<-bac_length))

figure(1)
hold on

[X,Y,Z] = cylinder(Rb,1000);
Z(2, :) = bac_length; Z(1, :) = -bac_length;

h = surf(X,Z,Y);
alpha 0.3
set(h,'edgecolor','none','facecolor','b')
hold on

[X,Y,Z] = sphere(1000);
X2 = X(:,1:500) * Rb;  X3 = X(:,500:end) * Rb;
Y2 = Y(:,1:500) * Rb;  Y3 = Y(:,500:end) * Rb;
Z2 = Z(:,1:500) * Rb;  Z3 = Z(:,500:end) * Rb;

h1 = surf(X2,Y2-bac_length,Z2);
alpha 0.3
set(h1,'edgecolor','none','facecolor','b')
h2 = surf(X3,Y3+bac_length,Z3);
alpha 0.3
set(h2,'edgecolor','none','facecolor','b')
hold on

            
plot3(MatP(p,1),MatP(p,2),MatP(p,3),'rd','MarkerSize',12,'MarkerFaceColor','r')
plot3(MatP(p,4),MatP(p,5),MatP(p,6),'gd','MarkerSize',12,'MarkerFaceColor','g')
%axis([-5 5 -5 5])
%axis([-4 4 -4 4 -4 4]*0.5)
%view(45, 45)
pause(0.1)

end

end