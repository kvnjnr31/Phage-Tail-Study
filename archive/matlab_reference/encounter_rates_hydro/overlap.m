%check for overlaps

function [] = overlap(MatP, MatB, D_para_y, D_para_z, Rb)

for p=1:size(MatP,1) 
r = sqrt((MatP(p,1)-MatB(1,1))^2 + D_para_y*(MatP(p,2)-MatB(1,2))^2 + D_para_z*(MatP(p,3)-MatB(1,3))^2);

%randomize locations
if r<(Rb)

figure(100)
hold on
plot(MatP(p,1),MatP(p,3),'ro','MarkerSize',12)
plot(MatP(p,7),MatP(p,9),'go','MarkerSize',12)
axis([-4 4 -4 4]*0.5)
%axis([-4 4 -4 4 -4 4]*0.5)
%view(45, 45)
pause(0.1)

end

end