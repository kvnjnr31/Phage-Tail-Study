for k = 1:1
    k
P = k*(1e+3);
B = 1;
nsteps = 1000;

%Parameters----------------------------
v=0; Rb=1; Rp=0.1; dt=0.01; Time=0; LX = 10; LY = 10;
x_left=-LX/2; x_right=LX/2; y_left=-LX/2; y_right=LX/2; V = LX*LX*LX; Rcrit=Rb+Rp;

omega = 1;                           %Number of Omega's
D0 = 2.54;                           %Diffusion Constant
nu = 0.82;                           %Diffusion Exponent
tmin = 0.001;                        %power-law gamma parameter
%***************************************************************************
%                            Main Loop
%***************************************************************************
kernel_new = zeros(1,nsteps);     %initialize kernel
en_indx = ones(1,P);   %initialize phage activation index.

%***************************************************************************
%                            Initialize System
%***************************************************************************
%initialize bacteria
for b=1:B

%fixed bacteria locations
MatB(b,1) = 0; MatB(b,2) = 0; %pos-new
MatB(b,3) = 0; MatB(b,4) = 0; %pos-old

MatB(1,5) = -v; MatB(1,6) = 0;

end

%initialize phage
for p=1:P 
    
r = 0;

while r<(Rb+Rp)
    
%random locations
x_old=(x_right-x_left-2*Rp)*rand + x_left + Rp;
y_old=(y_right-y_left-2*Rp)*rand + y_left + Rp;

%store the phage tracks
MatP(p,1) = x_old; MatP(p,2) = y_old;  %pos-new
MatP(p,3) = x_old; MatP(p,4) = y_old;  %pos-old

r = sqrt( (MatP(p,1)-MatB(1,1))^2 + (MatP(p,2)-MatB(1,2))^2);

end
%transition timer
MatP(p,5)=0;

%activate state
MatP(p,6)=1;

end

%***************************************************************************
%                     Simulate system over nsteps
%***************************************************************************
n = 1;
indx=0;

while n<nsteps

%Phage Dynamics
kk = find(MatP(:,6)==0);
jj = find(MatP(:,6)==1);

%Plotter-------------------------------------------------------------------
figure(1)
%
subplot(2,2,[1 3])
[Bx,By] = circle(MatB(:,1),MatB(:,2),Rb,100);
[Pxi,Pyi] = circle(MatP(kk,1),MatP(kk,2),Rp,100);
[Pxa,Pya] = circle(MatP(jj,1),MatP(jj,2),Rp,100);
plot(Bx,By,'k.',Pxi,Pyi,'r.',Pxa,Pya,'g.')
axis 'square'
axis([-LX/2 LX/2 -LY/2 LY/2])
drawnow

if n>1
subplot(2,2,[2 4])
plot(dt*[1:n-1],kernel_new(1:n-1),'k',dt*[1:n-1], pi*((Rb+Rp)^2)*v + 0*[1:n-1],'LineWidth',3)
axis 'square'
xlabel('Time: sec')
ylabel('Kernel: \beta')
axis square
%axis([0 10 0 20])
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
end

%--------------------------------------------------------------------------
%}    
    wx = zeros(P,1);
    wy = zeros(P,1);   
  
    kk = find (Time - 0*MatP(:,5)>=0);

    for p=1:length(kk) 
    r = 0;
    while r<(Rb+Rp)
    
    %random jump
    Ux = rand-0.5;
    wx = sqrt(24*D0*dt)*Ux;
    Uy = rand-0.5;
    wy = sqrt(24*D0*dt)*Uy;

    %check for overlap
    r = sqrt( (MatP(kk(p),1) + wx - MatB(1,1))^2 + (MatP(kk(p),2) + wy - MatB(1,2))^2);
    end
    MatP(kk(p),3) = MatP(kk(p),1); MatP(kk(p),4) = MatP(kk(p),2); 
    MatP(kk(p),1) = MatP(kk(p),1)+wx; MatP(kk(p),2) = MatP(kk(p),2)+wy; 
    end

    
    %generate waiting times
    U1 = rand(length(kk),1);
    tau = tmin*((1./U1).^(1./nu));
    MatP(kk,5) = Time + tau;                             %store waiting time
 
    %compute velocity field!!
    r = sqrt((MatP(:,1)-MatB(1,1)).^2+(MatP(:,2)-MatB(1,2)).^2);

    dot = MatB(1,5)*(MatP(:,1)-MatB(1,1)) + MatB(1,6)*(MatP(:,2)-MatB(1,2));  %dot product at each pt.(scalar).

    v_comp_x = (3*Rb./(4*r)).*(1+(Rb*Rb)./(3*r.*r))*MatB(1,5);
    v_comp_y = (3*Rb./(4*r)).*(1+(Rb*Rb)./(3*r.*r))*MatB(1,6);
    
    x_comp_x = (3*Rb./(4*r)).*(1./(r.*r) - (Rb*Rb)./(r.^4)).*(MatP(:,1)-MatB(1,1)).*dot;
    x_comp_y = (3*Rb./(4*r)).*(1./(r.*r) - (Rb*Rb)./(r.^4)).*(MatP(:,2)-MatB(1,2)).*dot;
    
    
    %Langevin x-dir
        MatP(:,3) = MatP(:,1);
        MatP(:,1) = MatP(:,3) - (v_comp_x+x_comp_x)*dt + MatB(1,5)*dt;

    %Langevin y-dir
        MatP(:,4) = MatP(:,2);
        MatP(:,2) = MatP(:,4) - (v_comp_y+x_comp_y)*dt + MatB(1,6)*dt;
        
%Compute Encounters--------------------------------------------------------
[MatP, MatB, Nenc, Nact] = encounter_fun(MatP, MatB, Rb+Rp);
kernel_new(n) = (Nenc/(Nact))*(V/dt);  %update kernel  

%apply B.C's---------------------------------------------------------------
[MatP] = boundarys(MatP, x_right, y_right);
    
%--------------------------------------------------------------------------
Time = Time + dt; %update time
n = n + 1; %update steps

end %time

if k == 1
    color = 'r';
elseif k == 5
    color = 'g';
elseif k == 9
    color = 'r';
end

%{
figure(1)
hold on
plot(dt*[1:15:nsteps],kernel_new(1:15:nsteps),'color',color,'LineWidth',3)
axis 'square'
xlabel('Time: sec')
ylabel('Kernel: \beta')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
axis([0 dt*nsteps 0 80])
drawnow
%}

end
