function sim_plotter(MatP,MatB,Rb,Rp,LX,LY,Time,N_dep)
%Rotate phage locations onto the x-y plane------------------------------------
MatR = zeros(size(MatP,1),3);                    %make space for rotations.
for p = 1:size(MatP,1)
v = MatP(p,1:3)';                                  %vector of phage points. 
theta = -atan(v(3)/v(2));                               %angle of rotation.
R = [1     0            0;                                %rotation matrix.
     0 cos(theta) -sin(theta);
     0 sin(theta)  cos(theta)];
vr = R*v;                                                  %rotated vector.
MatR(p,:) = vr';                                             %store vector.
end

%figure(100)
%hold on
%plot3(MatP(:,1),MatP(:,2),MatP(:,3),'k.','MarkerSize',10)

%Plot histogram of phage concentration distribution------------------------
r = sqrt(MatP(:,1).^2 + MatP(:,2).^2 + MatP(:,3).^2);
%jj = find(MatP(:,1)<0);
%r(jj)=[];

jj = find(r>2.5);
r(jj)=[];
[n,radius] = hist(r,40);

disc_size = (radius(2)-radius(1))/2;
volumes = (4/3)*pi*(radius+disc_size).^3 - (4/3)*pi*(radius-disc_size).^3;

figure(1)
hold on
plot(radius,n./(sum(n)*volumes),'r','LineWidth',2)
axis square
xlabel('Radius')
ylabel('Phage Concentration')
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow
%simplot-------------------------------------------------------------------
%
[xunit,yunit] = circle(0,0,Rb,1000);

figure(2)
hold on
plot3(MatR(:,1),MatR(:,2),MatR(:,3),'k.','MarkerSize',10)
plot(xunit, yunit,'m','LineWidth',2)
xlabel('x')
ylabel('y')
zlabel('z')
view([45,45])
axis([-LX/2 LX/2 -LY/2 LY/2 -LY/2 LY/2])
axis 'square'
%pause(0.1)
%}

%{
figure(1)
subplot(1,2,1)
[Bx,By] = circle(MatB(:,1),MatB(:,2),Rb,100);
[Pxi,Pyi] = circle(MatP(kk,1),MatP(kk,2),Rp,100);
[Pxa,Pya] = circle(MatP(jj,1),MatP(jj,2),Rp,100);
plot(Bx,By,'k.',Pxi,Pyi,'r.',Pxa,Pya,'g.')
axis 'square'
axis([-LX/2 LX/2 -LY/2 LY/2])
drawnow
%}

%Deposition plot-----------------------------------------------------------
%{
figure(1)
subplot(1,2,2)
hold on
plot(Time,N_dep,'r.','MarkerSize',10)
xlabel('Time')
ylabel('N')
axis 'square'
%pause(0.1)
%}
%{
%Contour plot--------------------------------------------------------------
figure(100)
subplot(1,2,1)
colormap(gray)
contourf(c{2}, c{1}, n,'LineColor','none');
axis 'square'
caxis([0 0.5])
colorbar('eastoutside')
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
%}

%scratch-------------------------------------------------------------------
%kk=find(MatP(:,5)==0);
%jj=find(MatP(:,5)==1);
%n = find(sqrt(MatP(:,1).^2 + MatP(:,2).^2+MatP(:,3).^2) <1.2);

%{
[n,c] = hist3(X, [10 10]);
 contour(c{2}, c{1}, n);

%smooth the contour!
newpoints = 100;
[xq,yq] = meshgrid(...
            linspace(min(min(xrow,[],2)),max(max(xrow,[],2)),newpoints ),...
            linspace(min(min(ycol,[],1)),max(max(ycol,[],1)),newpoints )...
          );
BDmatrixq = interp2(xrow,ycol,BDmatrix,xq,yq,'cubic');
[c,h]=contourf(xq,yq,BDmatrixq);
%}

%{
figure(1)
subplot(1,2,1)
[Bx,By] = circle(MatB(:,1),MatB(:,2),Rb,100);
[Pxi,Pyi] = circle(MatP(kk,1),MatP(kk,2),Rp,100);
[Pxa,Pya] = circle(MatP(jj,1),MatP(jj,2),Rp,100);
plot(Bx,By,'k.',Pxi,Pyi,'r.',Pxa,Pya,'g.')
axis 'square'
axis([-LX/2 LX/2 -LY/2 LY/2])
drawnow
%}
