%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            program to compute random walks on body surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MatP] = Langevin_fun(MatP, MatB, Rb, D0, k_bind, bac_length, dt)

indx1 = find(MatP(:,9)==1);                           %find attached phage.

indx2 = find(abs(MatP(indx1,2))<=bac_length);                        %body. 
indx3 = find(MatP(indx1,2)>bac_length);                              %head.
indx4 = find(MatP(indx1,2)<-bac_length);                             %tail.

%---------------------kick and project phage on surface--------------------
for kk=1:3
    
    if kk==1
        indx5 = indx1(indx2);           %attached phage within body limits.
        N = length(indx5);
    elseif kk==2
        indx5 = indx1(indx3);           %attached phage within head limits.
        N = length(indx5);
    elseif kk==3
        indx5 = indx1(indx4);           %attached phage within tail limits.
        N = length(indx5);
    end
    
%compute gradient at phage locations
r0 = [MatP(indx5,1), MatP(indx5,2), MatP(indx5,3)];

    if kk==1 
        n = [2*r0(:,1), zeros(N,1), 2*r0(:,3)];             %body gradient.
    elseif kk==2
        n = [2*r0(:,1), 2*(r0(:,2)-bac_length), 2*r0(:,3)]; %head gradient.  
    elseif kk==3
        n = [2*r0(:,1), 2*(r0(:,2)+bac_length), 2*r0(:,3)]; %tail gradient. 
    end
    
%-----------------------------random kick----------------------------------
w = rand(N,3);                                      %uniform random vector.
%random vec orthogonal to both n and w.
v = [n(:,2).*w(:,3)-n(:,3).*w(:,2), -(n(:,1).*w(:,3)-n(:,3).*w(:,1)), n(:,1).*w(:,2)-n(:,2).*w(:,1)];

Ux = rand(N,1)-0.5; L = sqrt(24*D0*dt).*Ux;
v = L.*v./(sqrt(sum(v.*v,2)));           %kick in direction of rand vector.

%rotate random kick
nxv = [n(:,2).*v(:,3)-n(:,3).*v(:,2), -(n(:,1).*v(:,3)-n(:,3).*v(:,1)), n(:,1).*v(:,2)-n(:,2).*v(:,1)];
ndotv = sum(n.*v,2);
theta = 2*pi*rand(N,1);                        %angle between zero and 2pi.
%rotated random vector.
v = v.*cos(theta) + nxv.*sin(theta) + n.*(ndotv).*(1-cos(theta)); 

%----------------------------Force Vectors---------------------------------
F_proj = zeros(N,3);
i_bind = MatP(indx5,8);                           %phage binding receptors.

r_bind = MatB(i_bind(i_bind>0),1:3);     %location of binding sites (i_bind>0).

v_bind = r_bind - r0(i_bind>0,:);       %binding vector (points to receptor).
v_bind_mag = sqrt(sum(v_bind.*v_bind,2));          %length of binding vecs.

%force lengths.
x = v_bind_mag - MatP(indx5(i_bind>0),10) - MatP(indx5(i_bind>0),11);
H = (1 + x/abs(x))/2;                                  %heavyside function.
F = H.*k_bind*x.*v_bind./v_bind_mag;                          %Binding forces.
F_proj(i_bind>0,:) = F - sum(F.*n(i_bind>0,:),2).*n(i_bind>0,:)./sum(n(i_bind>0,:).*n(i_bind>0,:),2);     %force vectors projected onto plane.

%kicked vector with Forces!
r1 = r0 + v + 0*F_proj;  

%project kicked vector back towards surface
phi_r1 = zeros(N,1);
grad_r1 = zeros(N,3);

indx6 = find(abs(r1(:,2))<=bac_length);                              %body. 
indx7 = find(r1(:,2)>bac_length);                                    %head.
indx8 = find(r1(:,2)<-bac_length);                                   %tail.

for jj=1:3
    
    if jj==1
        tail_length = MatP(indx5(indx6),10);
        Rp = MatP(indx5(indx6),11);
        phi_r1(indx6) = r1(indx6,1).^2 + zeros(length(indx6),1) + r1(indx6,3).^2 - (Rb+(tail_length+Rp)/2).^2;
        grad_r1(indx6,:) = [2*r1(indx6,1), zeros(length(indx6),1), 2*r1(indx6,3)]; 
    elseif jj==2
        tail_length = MatP(indx5(indx7),10);
        Rp = MatP(indx5(indx7),11);
        phi_r1(indx7) = r1(indx7,1).^2 + (r1(indx7,2)-bac_length).^2  + r1(indx7,3).^2 - (Rb+(tail_length+Rp)/2).^2;
        grad_r1(indx7,:) = [2*r1(indx7,1), 2*(r1(indx7,2)-bac_length), 2*r1(indx7,3)];
    elseif jj==3
        tail_length = MatP(indx5(indx8),10);
        Rp = MatP(indx5(indx8),11);
        phi_r1(indx8) = r1(indx8,1).^2 + (r1(indx8,2)+bac_length).^2 + r1(indx8,3).^2 - (Rb+(tail_length+Rp)/2).^2;
        grad_r1(indx8,:) = [2*r1(indx8,1), 2*(r1(indx8,2)+bac_length), 2*r1(indx8,3)]; 
    end
    
end

r2 = r1 - phi_r1.*grad_r1./(sum(grad_r1.*grad_r1,2)); %projected vec!
  
%Langevin x-dir
MatP(indx5,4) = 0*MatP(indx5,1);                           %x-drift velocity.
MatP(indx5,1) = r2(:,1); 
                                      
%Langevin y-dir
MatP(indx5,5) = 0*MatP(indx5,2);                  %y-drift velocity.
MatP(indx5,2) = r2(:,2);
                                     
%Langevin z-dir
MatP(indx5,6) = 0*MatP(indx5,3);                  %z-drift velocity.
MatP(indx5,3) = r2(:,3);

end
