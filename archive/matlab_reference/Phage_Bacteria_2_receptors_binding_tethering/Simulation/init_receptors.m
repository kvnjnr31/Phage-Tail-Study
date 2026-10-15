%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Algorithm to initialize bcteria receptors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MatB] = init_receptors(B, delta, m, Rb, bac_length, Rr)
%Randomize locations 
theta_max = 2*(bac_length+Rb-(delta/cos(m)))/m;

for kk = 1:B    
output = 1; 

    while output == 1    
    a = 2*delta*rand - delta; theta = theta_max*rand;
    y = m*theta + (a+delta)/cos(m) - bac_length - Rb ;
    
        if abs(y)<bac_length    %body receptor
            MatB(kk,1) = Rb*cos(theta); MatB(kk,2) = y; MatB(kk,3) = Rb*sin(theta);
        elseif y>bac_length   %head receptor
            MatB(kk,2) = y;
            y = y - bac_length;
            MatB(kk,1) = sqrt(Rb^2-y*y).*cos(theta); MatB(kk,3) = sqrt(Rb^2-y*y).*sin(theta);
        elseif y<bac_length   %tail receptor
            MatB(kk,2) = y;
            y = y + bac_length;
            MatB(kk,1) = sqrt(Rb^2-y*y).*cos(theta); MatB(kk,3) = sqrt(Rb^2-y*y).*sin(theta);
        end
        
    [output] = receptor_overlap(MatB, kk, Rr);

    end   
end
