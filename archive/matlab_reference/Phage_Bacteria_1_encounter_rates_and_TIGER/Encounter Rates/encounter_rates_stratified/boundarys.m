function [MatP, MatB] = boundarys(MatP, MatB, x_right, y_right, z_right)
%Apply random boundary conditions to phage

    jj = find(abs(MatP(:,1))>=x_right);
    %
    MatP(jj,1) = MatP(jj,1) - (2*(x_right))*MatP(jj,1)./abs(MatP(jj,1));                                      
    MatP(jj,7) = MatP(jj,1);
    
    MatP(jj,5) = 1;                                              %activate.
    %}
        
    jj = find(abs(MatP(:,2))>=y_right);
    %
    MatP(jj,2) = MatP(jj,2) - (2*(y_right))*MatP(jj,2)./abs(MatP(jj,2)); 
    MatP(jj,8) = MatP(jj,2);
    
    
    MatP(jj,5) = 1;                                              %activate.
    %}
    
    jj = find(abs(MatP(:,3))>=z_right);
    %
    MatP(jj,3) = MatP(jj,9);
    
    MatP(jj,5) = 1;  
    
%Apply random boundary conditions to bacteria

    jj = find(abs(MatB(:,1))>=x_right);
    %
    MatB(jj,1) = MatB(jj,1) - (2*(x_right))*MatB(jj,1)./abs(MatB(jj,1));                                      
    MatB(jj,7) = MatB(jj,1);
    
    jj = find(abs(MatB(:,2))>=y_right);
    %
                                         
    MatB(jj,2) = MatB(jj,2) - (2*(y_right))*MatB(jj,2)./abs(MatB(jj,2)); 
    MatB(jj,8) = MatB(jj,2);
    
    
    jj = find(abs(MatB(:,3))>=z_right);
    %
    MatB(jj,1) = MatB(jj,7);
    MatB(jj,2) = MatB(jj,8);
    MatB(jj,3) = MatB(jj,9);
    
    %MatB(jj,4) = -MatB(jj,4);
    %MatB(jj,5) = -MatB(jj,5);
    MatB(jj,6) = -MatB(jj,6);
    
end
