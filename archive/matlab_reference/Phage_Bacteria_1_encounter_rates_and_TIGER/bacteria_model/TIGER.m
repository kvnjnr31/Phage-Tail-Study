% Algorithm to simulate swimming bacteria
% parameters -----------------------------
h = 0.77;
kE = 1.3;
k = 1.3;
af = 0.1;
w = -100;
Rb = 1;

%Movie Stuff---------------------------------------------------------------
%tfinal=Sim_Length*dt;%500;
%Nsnaps=1;
imovie=0;
savemovie=1;
clear M
%kstep=(tfinal/(1*Nsnaps));

%flagellum ---------------------------------
x = linspace(0,7,25);
Ex = 1 - exp(-(kE^2)*x.^2);

t=0;
dt = 0.05;

% Generate initial traction points;
head_coord_x0 =  linspace(-2*(3*Rb),0,10);
head_coord_x = head_coord_x0;
hy = 0*head_coord_x0;
hz =  (1*Rb)*sqrt(1-((head_coord_x0+(3*Rb))/((3*Rb))).^2);
head_theta0 = [hy;hz];
head_coord_yz = head_theta0;
theta = pi/3; %rotation angle
R_theta = [cos(theta) -sin(theta); sin(theta) cos(theta)];

for ii=1:6
head_theta0 = R_theta*head_theta0;
head_coord_yz = [head_coord_yz head_theta0];
head_coord_x = [head_coord_x head_coord_x0];
end

theta = -pi/24; %rotation angle
R_theta = [cos(theta) -sin(theta); sin(theta) cos(theta)];

for jj = 1:100
 
   head_coord_yz = R_theta*head_coord_yz;    
    
    y = h*Ex.*(sin(k*x-w*t));
    z = h*Ex.*(cos(k*x-w*t));
  
    t = t + dt;
 
 %---------------------------------
 imovie=imovie+1; 
 %---------------------------------
 
 figure(1)
 plot3(head_coord_x,head_coord_yz(1,:),head_coord_yz(2,:),'ko--','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','k')
 hold on
 plot3(x,y,z,'ko--','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','k')
 
 grid off;
 title(['TIGER Time: ' num2str(t,2) ' (sec)'])
 xlabel('L_x')
 ylabel('L_y')
 zlabel('L_z')
 set(gca,'xlim',[-2*(3*Rb) 7],'ylim',[-2 2],'zlim',[-2 2])
 view([10,-22,2])
 axis 'square'
 set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
 drawnow
 hold off
 
 %savemovie-----------------------------------------------------------------
 g=figure(1);
 if(savemovie)
 M(imovie)=getframe(g);
 end
    
end

%savemovie----------------------------------------------------------------
movie2avi(M, 'TIGER.avi', 'compression', 'None');



