opengl software           %to speed up plotting

for k = 5:5
    
k
P = k*(1e+4);
B = 1;
nsteps = 100;

%Parameters----------------------------
v=30; Rb=1; Rp=0.1; dt=0.01; Time=0; LX = 50; LY = 10;
x_left=-LX/2; x_right=LX/2; y_left=-LX/2; y_right=LX/2; z_left=-LX/2; z_right=LX/2; 
V = LX*LX*LX; Rcrit=Rb+Rp;

omega = 1;                                              %Number of Omega's.
D0 = 2.54;                                             %Diffusion Constant.
kernel_new = zeros(1,nsteps);                           %initialize kernel.

%**************************************************************************
%                            Initialize System
%**************************************************************************

%initialize bacteria
for b=1:B
%fixed bacteria locations
MatB(b,1) = 0; MatB(b,2) = 0; MatB(b,3) = 0;                      %pos-new.
MatB(1,4) = -v;  MatB(1,5) = 0;   MatB(1,6) = 0;                 %velocity.
end

%initialize phage
for p=1:P 
r = 0;

    while r<(Rb+Rp)         %ensure no particles are within the capture radius.
    
        %random locations
        x_old=(x_right-x_left-2*Rp)*rand + x_left + Rp;
        y_old=(y_right-y_left-2*Rp)*rand + y_left + Rp;
        z_old=(z_right-z_left-2*Rp)*rand + z_left + Rp;

        %store the phage tracks
        MatP(p,1) = x_old; MatP(p,2) = y_old;  MatP(p,3) = z_old; %pos-new
        r = sqrt( (MatP(p,1)-MatB(1,1))^2 + (MatP(p,2)-MatB(1,2))^2 + (MatP(p,3)-MatB(1,3))^2);

    end

%transition timer
MatP(p,4)=0;

%activate state
MatP(p,5)=1;

end

%**************************************************************************
%                            Main Loop
%**************************************************************************
n = 1;
indx=0;

while n<nsteps

%Phage Dynamics
%
kk = find(MatP(:,5)==0);
jj = find(MatP(:,5)==1);

%Plotter-------------------------------------------------------------------
figure(1)
%{
subplot(2,2,[1 3])
[Bx,By] = circle(MatB(:,1),MatB(:,2),Rb,100);
[Pxi,Pyi] = circle(MatP(kk,1),MatP(kk,2),Rp,100);
[Pxa,Pya] = circle(MatP(jj,1),MatP(jj,2),Rp,100);
plot(Bx,By,'k.',Pxi,Pyi,'r.',Pxa,Pya,'g.')
axis 'square'
axis([-LX/2 LX/2 -LY/2 LY/2])
drawnow
%}
if n>1
subplot(2,2,[2 4])
plot(dt*[1:n-1],kernel_new(1:n-1),'k',dt*[1:n-1], pi*((Rb+Rp)^2)*v + 0*[1:n-1],'LineWidth',3)
axis 'square'
xlabel('Time: sec')
ylabel('Kernel: \beta')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
end

%{
figure(2)
plot3(MatP(:,1),MatP(:,2),MatP(:,3),'ko' )
xlabel('x')
ylabel('y')
ylabel('z')
%}
%--------------------------------------------------------------------------
%}    
   
%Langevin x-dir
    MatP(:,1) = MatP(:,1) + MatB(1,4)*dt;

%Langevin y-dir
    MatP(:,2) = MatP(:,2) + MatB(1,5)*dt;
    
%Langevin z-dir
    MatP(:,3) = MatP(:,3) +  MatB(1,6)*dt;
        
%Compute Encounters--------------------------------------------------------
[MatP, MatB, Nenc, Nact] = encounter_fun_3D(MatP, MatB, Rcrit);
kernel_new(n) = (Nenc/(Nact))*(V/dt);  %update kernel  

%apply B.C's---------------------------------------------------------------
[MatP] = boundarys_3D(MatP, x_right, y_right, z_right, Rp);
    
%--------------------------------------------------------------------------
Time = Time + dt; %update time
n = n + 1; %update steps

end %time

end
