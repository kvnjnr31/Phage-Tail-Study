%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Program to integrate Langevin Eqns
%%Input parameters
%MatP MatB: Phage, bacteria matrix
%dt: diffusion time
%wx, wy, wz: noise

%Output parameters
%MatP MatB: Phage, bacteria matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MatP, MatB] = Langevin_fun_3D(MatP, MatB, dt, wx, wy, wz)

    %Langevin x-dir
    MatP(:,7) = MatP(:,1);                               %x-drift velocity.
    MatP(:,1) = MatP(:,1) + wx;
    
    MatB(:,8) = MatB(:,1);
    MatB(:,1) = MatB(:,1) + MatB(:,4)*dt;
    
    %Langevin y-dir
    MatP(:,8) = MatP(:,2);                               %y-drift velocity.
    MatP(:,2) = MatP(:,2) + wy;
    
    MatB(:,9) = MatB(:,2);
    MatB(:,2) = MatB(:,2) + MatB(:,5)*dt;
    
    %Langevin z-dir
    MatP(:,9) = MatP(:,3);                      %z-drift velocity.
    MatP(:,3) = MatP(:,3) + wz;

    MatB(:,10) = MatB(:,3);
    MatB(:,3) = MatB(:,3) + MatB(:,6)*dt;