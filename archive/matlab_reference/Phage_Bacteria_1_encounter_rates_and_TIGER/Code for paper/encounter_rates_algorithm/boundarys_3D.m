%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Program to apply random boundary conditions to phage and bacteria
%Input parameters
%MatP MatB: Phage, bacteria matrix
%x_right, y_right, z_right: boundary locations 
%Rp, Rb: phage and bacteria radii

%Output parameters
%MatP MatB: Phage, bacteria matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MatP,MatB] = boundarys_3D(MatP, MatB, x_right, y_right, z_right, Rp, Rb)

%phage---------------------------------------------------------------------

    jj = find(abs(MatP(:,1))>=x_right);
    MatP(jj,1) = MatP(jj,1) - (2*(x_right))*MatP(jj,1)./abs(MatP(jj,1));
    MatP(jj,2) = 2*(y_right-Rp)*rand(length(jj),1) - y_right;  %rand y-dir.
    MatP(jj,3) = 2*(z_right-Rp)*rand(length(jj),1) - z_right;  %rand z-dir.
    
    jj = find(abs(MatP(:,2))>=y_right);
    MatP(jj,1) = 2*(y_right-Rp)*rand(length(jj),1) - y_right;  %rand y-dir.
    MatP(jj,2) = MatP(jj,2) - (2*(y_right))*MatP(jj,2)./abs(MatP(jj,2));
    MatP(jj,3) = 2*(z_right-Rp)*rand(length(jj),1) - z_right;  %rand z-dir.
    
    jj = find(abs(MatP(:,3))>=z_right);
    MatP(jj,1) = 2*(y_right-Rp)*rand(length(jj),1) - y_right;  %rand y-dir.
    MatP(jj,2) = 2*(z_right-Rp)*rand(length(jj),1) - z_right;  %rand z-dir.
    MatP(jj,3) = MatP(jj,3) - (2*(z_right))*MatP(jj,3)./abs(MatP(jj,3));

%bacteria------------------------------------------------------------------

    jj = find(abs(MatB(:,1))>=x_right);
    MatB(jj,1) = MatB(jj,1) - (2*(x_right))*MatB(jj,1)./abs(MatB(jj,1));
    MatB(jj,2) = 2*(y_right-Rb)*rand(length(jj),1) - y_right;  %rand y-dir.
    MatB(jj,3) = 2*(z_right-Rb)*rand(length(jj),1) - z_right;  %rand z-dir.
    
    jj = find(abs(MatB(:,2))>=y_right);
    MatB(jj,1) = 2*(y_right-Rb)*rand(length(jj),1) - y_right;  %rand y-dir.
    MatB(jj,2) = MatB(jj,2) - (2*(y_right))*MatB(jj,2)./abs(MatB(jj,2));
    MatB(jj,3) = 2*(z_right-Rb)*rand(length(jj),1) - z_right;  %rand z-dir.
    
    jj = find(abs(MatB(:,3))>=z_right);
    MatB(jj,1) = 2*(y_right-Rb)*rand(length(jj),1) - y_right;  %rand y-dir.
    MatB(jj,2) = 2*(z_right-Rb)*rand(length(jj),1) - z_right;  %rand z-dir.
    MatB(jj,3) = MatB(jj,3) - (2*(z_right))*MatB(jj,3)./abs(MatB(jj,3));
                                      
end
