function [] = RATTLE(MatP, MatB, bac_length, dt)

N = size(MatP,1); M = size(MatB,1);

%Lennard Jones parameters
kBT=300*(1.38e-23)*(1e12);                            %thermal temperature.
sigma = 1e-4; epsilon = 0.5*kBT;   
rc = 1;             %r critical
mass = (0.1)*1e-15; %using a density of water

for ii = 1:N
    
    xm = [MatP(ii,1); MatP(ii,2); MatP(ii,3)]; %phage location: column vec!
    xm_plus = xm;
    pm = [MatP(ii,4); MatP(ii,5); MatP(ii,6)];       %momentum: column vec!
    
    lambda = 0; r_mag = 2; eta = 1e-4;
    
    %compute normals to the surface --------------------------------------- 
    if(abs(xm(2))<=bac_length) %body. 
        nm = [2*xm(1); 0; 2*xm(3)];                       
    elseif(xm(2)>bac_length) %head.
        nm = [2*xm(1); 2*(xm(2)-bac_length); 2*xm(3)];    
    elseif(xm(2)<-bac_length) %tail.
        nm = [2*xm(1); 2*(xm(2)+bac_length); 2*xm(3)];  
    end
    
    %compute total forces-------------------------------------------------- 
    d2 = ((xm(1)-MatB(:,1)).^2 + (xm(2)-MatB(:,2)).^2 + ...
          (xm(3)*MatB(:,3)).^2);           %vector distance from receptors.
     
    dr_dx = 1./sqrt(d2).*[(xm(1)-MatB(:,1)), ...
          (xm(2)-MatB(:,2)), (xm(3)-MatB(:,3))];
        
    dLJ_dx = 4*epsilon*[12*(sigma./sqrt(d2)).^11 - ...
             6*(sigma./sqrt(d2)).^5].*(-sigma./(sqrt(d2).^2)).*dr_dx;
    
    jj = find(sqrt(d2)<rc);
    
    f = sum(dLJ_dx(jj,:),1);                  %total force in each channel.
    f = f';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                              Main Loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    while r_mag < eta     
        %compute normals to the surface------------------------------------ 
        if(abs(xm_plus(2))<=bac_length) %body.                     
            nm_plus = [2*xm_plus(1); 0; 2*xm_plus(3)]; 
        elseif(xm_plus(2)>bac_length) %head. 
            nm_plus = [2*xm_plus(1); 2*(xm_plus(2)-bac_length); 2*xm_plus(3)];  
        elseif(xm_plus(2)<-bac_length) %tail.
            nm_plus = [2*xm_plus(1); 2*(xm_plus(2)+bac_length); 2*xm_plus(3)]; 
        end

        %compute momentum--------------------------------------------------
        pm_plus_half = pm + dt/2*(f - lambda*nm);

        %compute constraint------------------------------------------------
        if(abs(xm_plus(2))<=bac_length) %body. 
            gm_plus = xm_plus(1)^2 + 0 + xm_plus(3)^2 - Rb^2;                
        elseif(xm_plus(2)>bac_length) %head.
            gm_plus = xm_plus(1)^2 + (xm_plus(2)-bac_length).^2 + xm_plus(3)^2 - Rb^2; 
        elseif(xm_plus(2)<-bac_length) %tail.
            gm_plus = xm_plus(1)^2 + (xm_plus(2)+bac_length).^2 + xm_plus(3)^2 - Rb^2; 
        end

        %compute r---------------------------------------------------------
        r = [xm - xm_plus + (dt/mass)*pm_plus_half; gm_plus];

        %compute Jacobi----------------------------------------------------
        J = [-1      0       0      -(dt^2/(2*mass))*nm(1);
             0      -1       0      -(dt^2/(2*mass))*nm(2);
             0      0       -1      -(dt^2/(2*mass))*nm(3);
             nm(1)  nm(2)   nm(3)   0];

        %update position---------------------------------------------------
        x_lambda = [xm_plus;lambda] - inv(J)*r;
        
        %compute new r_mag
        lambda = x_lambda(4);
        r_mag = sqrt(sum(r.*r));
    end  
    
    xm_plus = x_lambda(1:3);
    MatP(ii,1) = xm_plus(1); MatP(ii,2) = xm_plus(2); MatP(ii,3) = xm_plus(3);
    MatP(ii,4) = pm_plus_half(1); MatP(ii,5) = pm_plus_half(2); MatP(ii,6) = pm_plus_half(3);       
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Update Velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:N

    xm_plus = [MatP(ii,1); MatP(ii,2); MatP(ii,3)];
            
    %compute total forces-------------------------------------------------- 
    d2 = ((xm_plus(1)-MatB(:,1)).^2 + (xm_plus(2)-MatB(:,2)).^2 + ...
          (xm_plus(3)*MatB(:,3)).^2);      %vector distance from receptors.
     
    dr_dx = 1./sqrt(d2).*[(xm_plus(1)-MatB(:,1)), ...
          (xm_plus(2)-MatB(:,2)), (xm_plus(3)-MatB(:,3))];
        
    dLJ_dx = 4*epsilon*[12*(sigma./sqrt(d2)).^11 - ...
             6*(sigma./sqrt(d2)).^5].*(-sigma./(sqrt(d2).^2)).*dr_dx;
    
    jj = find(sqrt(d2)<rc);
    
    f = sum(dLJ_dx(jj,:),1);                  %total force in each channel.
    fm_plus = f';
    
    %compute normals to the surface------------------------------------ 
    if(abs(xm_plus(2))<=bac_length) %body.                     
        nm_plus = [2*xm_plus(1); 0; 2*xm_plus(3)]; 
    elseif(xm_plus(2)>bac_length) %head. 
        nm_plus = [2*xm_plus(1); 2*(xm_plus(2)-bac_length); 2*xm_plus(3)];  
    elseif(xm_plus(2)<-bac_length) %tail.
        nm_plus = [2*xm_plus(1); 2*(xm_plus(2)+bac_length); 2*xm_plus(3)]; 
    end
    
    %update momentum-------------------------------------------------------
    pm_plus_half = [MatP(ii,4); MatP(ii,5); MatP(ii,6)];       %momentum: column vec!
    pm_plus = pm_plus_half + (dt/2)*fm_plus; 
    mu = 0; r_mag = 2; eta = 1;
    
    while r_mag < eta
     
        %compute acceleration
        dpdt_m_plus = fm_plus - mu*nm_plus;
        
        %compute r
        r = [pm_plus_half - pm_plus + (dt/2)*dpdt_m_plus; sum(pm_plus.*nm_plus)/mass];
        
        %compute Jacobi----------------------------------------------------
        J = [-1             0           0       -(dt^2/(2))*nm_plus(1);
             0              -1          0       -(dt^2/(2))*nm_plus(2);
             0              0           -1      -(dt^2/(2))*nm_plus(3);
             nm(1)/mass  nm(2)/mass   nm(3)/mass   0];
         
        %update position---------------------------------------------------
        p_mu = [pm_plus;mu] - inv(J)*r;
        
        %compute new r_mag
        mu = p_mu(4);
        r_mag = sqrt(sum(r.*r));
        
    end
end


    
    
    
    
    
    
    
    

