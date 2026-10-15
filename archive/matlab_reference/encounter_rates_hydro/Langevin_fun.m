function [MatP] = Langevin_fun(MatP, MatB, Rb, dt, D_para_y, D_para_z)

%-----------------program to compute Stokes flow around--------------------
%compute velocity field!!
    r = sqrt((MatP(:,1)).^2+D_para_y*(MatP(:,2)).^2 + D_para_z*(MatP(:,3)).^2);
    r2 = (r.*r);
    r3 = (r.*r.*r);
    r4 = (r.*r.*r.*r);
    
    %Stokes Flow
    dot = MatB(1,4)*(MatP(:,1)) + D_para_y*MatB(1,5)*(MatP(:,2)) + D_para_z*MatB(1,6)*(MatP(:,3));  %dot product at each pt.(scalar).

    v_comp_x = (3*Rb./(4*r)).*(1+(Rb*Rb)./(3*r2))*MatB(1,4);
    v_comp_y = D_para_y*(3*Rb./(4*r)).*(1+(Rb*Rb)./(3*r2))*MatB(1,5);
    v_comp_z = D_para_z*(3*Rb./(4*r)).*(1+(Rb*Rb)./(3*r2))*MatB(1,6);
    
    x_comp_x = (3*Rb./(4*r)).*(1./(r2) - (Rb*Rb)./r4).*(MatP(:,1)).*dot;
    x_comp_y = D_para_y*(3*Rb./(4*r)).*(1./(r2) - (Rb*Rb)./r4).*(MatP(:,2)).*dot;
    x_comp_z = D_para_z*(3*Rb./(4*r)).*(1./(r2) - (Rb*Rb)./r4).*(MatP(:,3)).*dot;
    
    %Rotational Flow
    phi_2 = -1*(Rb*Rb*Rb)*MatP(:,2)./r3;
    phi_3 = -1*(Rb*Rb*Rb)*MatP(:,3)./r3;
    
    %Langevin x-dir
    MatP(:,7) = MatP(:,1);                               %x-drift velocity.
    MatP(:,1) = MatP(:,1) - (v_comp_x+x_comp_x)*dt + MatB(1,4)*dt;
                                      
    %Langevin y-dir
    MatP(:,8) = D_para_y*MatP(:,2);                      %y-drift velocity.
    MatP(:,2) = D_para_y*(MatP(:,2) - (v_comp_y+x_comp_y)*dt - MatB(1,4)*phi_3*dt + MatB(1,5)*dt);
                                     
    %Langevin z-dir
    MatP(:,9) = D_para_z*MatP(:,3);                      %z-drift velocity.
    MatP(:,3) = D_para_z*(MatP(:,3) - (v_comp_z+x_comp_z)*dt + MatB(1,4)*phi_2*dt + MatB(1,6)*dt);
