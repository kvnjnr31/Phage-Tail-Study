%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Algorithm to initialize phage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [MatP] = init_phage(P, Rb, Rp, x_right, x_left, ...
          y_right, y_left, z_right, z_left, Ltail)

for p=1:P 
    
MatP(p,10)= Ltail;                                            %tail length.
MatP(p,11) = Rp;                                              %head radius.

%randomize locations-------------------------------------------------------
r = 0;
    while (r<(Rb+Rp))||(r>(Rb+Rp+0.05))
    x_old=  ((x_right-x_left-2*Rp)*rand + x_left + Rp); 
    y_old=  ((y_right-y_left-2*Rp)*rand + y_left + Rp); 
    z_old=  ((z_right-z_left-2*Rp)*rand + z_left + Rp);
    r = sqrt((x_old)^2 + (z_old)^2);
    end
%--------------------------------------------------------------------------
%store the phage tracks
MatP(p,1) = x_old; MatP(p,2) = y_old;  MatP(p,3) = z_old; 

%initial xyz-dir positions 
MatP(p,4)=x_old;  MatP(p,5)=y_old; MatP(p,6)=z_old;

%transition timer
MatP(p,7)=0;
%activate state
MatP(p,8)=0;
%surface attachment indicator
MatP(p,9)=0;

end
