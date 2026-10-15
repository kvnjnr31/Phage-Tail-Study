%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Program to find jumping phage and generate Brownian noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MatP] = noise_gen(MatP, D0, dt, tmin, nu, Time, Rb, bac_length)

dd = find(MatP(:,9) == 0);                         %find detached phage.
jj = find (Time - MatP(dd,7)>=0);       %find detached phage ready to jump.

dj = dd(jj);                                  %indices of detached jumpers.

%generate waiting times
U1 = rand(length(dj),1); tau = tmin*((1./U1).^(1./nu));
MatP(dj,7) = Time + tau;                               %store waiting time
     
   
for m = 1:length(dj)
    
    tail_length = MatP(dj(m),10);
    Rp = MatP(dj(m),11);
    
    phi = 0;
    phi1 = 0;
    phi2 = 0;

    n = 0;

    while (phi<=0)||(phi1<=0)||(phi2<=0)              
        Ux = rand-0.5; wx = sqrt(24*D0*dt)*Ux;
        Uy = rand-0.5; wy = sqrt(24*D0*dt)*Uy;
        Uz = rand-0.5; wz = sqrt(24*D0*dt)*Uz;

        rx_new = MatP(dj(m),1) +  wx;
        ry_new = MatP(dj(m),2) +  wy;
        rz_new = MatP(dj(m),3) +  wz;
        
        if (abs(ry_new)<=bac_length)
        phi = rx_new^2 + rz_new^2 - (Rb+Rp)^2;
        else
          phi = 1;
        end
        phi1 = rx_new^2 + (ry_new-bac_length)^2 + rz_new^2 - (Rb+Rp)^2;
        phi2 = rx_new^2 + (ry_new+bac_length)^2 + rz_new^2 - (Rb+Rp)^2;

        n=n+1;
        %
        if n>1e+5
            figure(200)
            hold on
            
            if n==(1+ 1e+5)
            [X,Y,Z] = cylinder(Rb+Rp,1000);
            Z(2, :) = bac_length; Z(1, :) = -bac_length;

            h = surf(X,Z,Y);
            alpha 0.3
            set(h,'edgecolor','none','facecolor','b')
            hold on

            [X,Y,Z] = sphere(1000);
            X2 = X(:,1:500) * (Rb+Rp);  X3 = X(:,500:end) * (Rb+Rp);
            Y2 = Y(:,1:500) * (Rb+Rp);  Y3 = Y(:,500:end) * (Rb+Rp);
            Z2 = Z(:,1:500) * (Rb+Rp);  Z3 = Z(:,500:end) * (Rb+Rp);

            h1 = surf(X2,Y2-bac_length,Z2);
            alpha 0.3
            set(h1,'edgecolor','none','facecolor','b')
            h2 = surf(X3,Y3+bac_length,Z3);
            alpha 0.3
            set(h2,'edgecolor','none','facecolor','b')
            hold on
            end
            
            plot3(rx_new,ry_new,rz_new,'ko','MarkerSize',12,'MarkerFaceColor','k')
            plot3(MatP(dj(m),1),MatP(dj(m),2),MatP(dj(m),3),'rd','MarkerSize',12,'MarkerFaceColor','r')
            plot3(MatP(dj(m),4),MatP(dj(m),5),MatP(dj(m),6),'gd','MarkerSize',12,'MarkerFaceColor','g')
            xlabel('x')
            ylabel('y')
            zlabel('z')
            drawnow
        end
        %}
    end
    MatP(dj(m),4) = MatP(dj(m),1);
    MatP(dj(m),1) = rx_new;
    MatP(dj(m),5) = MatP(dj(m),2);
    MatP(dj(m),2) = ry_new;
    MatP(dj(m),6) = MatP(dj(m),3);
    MatP(dj(m),3) = rz_new;

    if ((n>1))
       MatP(dj(m),9)=1;    %attach phage 
    end           
end
 
    
       