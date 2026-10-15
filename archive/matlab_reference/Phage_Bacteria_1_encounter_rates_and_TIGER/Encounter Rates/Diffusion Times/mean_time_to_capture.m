opengl software

%Mean time to capture model
delta = 1;                                                    %jump length.
N = 100000;                                                  %number of steps.
P = 1e+3;                                              %number of particles.
MatP = zeros(P,4); 
%set activation
MatP(:,2) = ones(P,1);
Time = 0; dt = 0.01; 
LX = 10; %boundary

%Diffusion Parameters
nu=0.82;
D0 = (delta^2)/(2*dt);
tmin=dt*sinc(nu)^(1/nu);
Nact = 1; n=1;

%generate waiting times
U1 = rand(P,1);
tau = tmin*((1./U1).^(1./nu));
MatP(:,4) = tau;  

while(Nact>0)

%find jumpers and plot histogram
kk = find(Time - MatP(:,4)>=0);
length(kk)/P
Prob = nu*(tmin^(nu+1))/(dt^nu)

[nout nbins]=hist(MatP(:,4),50);
figure(99)
hold on
plot(Time,length(kk)/P,'k.',Time,Prob,'r.','LineWidth',2)
drawnow


Ux = rand(length(kk),1) - 0.5;                %uniform between -0.5 and 0.5.
wx = sqrt(24*D0*dt)*Ux;

MatP(kk,1) = MatP(kk,1) + wx;    %step to the left or right.

%generate waiting times
U1 = rand(length(kk),1);
tau = tmin*((1./U1).^(1./nu));
MatP(kk,4) = Time + tau;                             %store waiting time
    
%check for crossing boundaries
k = find(abs(MatP(:,1))>=LX);
j = find(MatP(k,2)>0);
MatP(k(j),3)=Time; %adsorption time.
MatP(k(j),2)=0; %deactivate.

%{
if(mod(n,100)==0)
figure(1)
hold on
plot(Time,mean(MatP(:,3)),'ko',[0:0.01:Time], (LX^2)/(2*nu*D0*(1-(tmin/(dt+tmin))^(nu))) + 0*[0:0.01:Time],'k')
drawnow
end
%}

Nact = sum(MatP(:,2));
Time = Time + dt;
n = n + 1;
end





