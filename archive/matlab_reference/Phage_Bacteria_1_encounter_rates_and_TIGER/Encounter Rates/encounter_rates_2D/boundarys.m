function [MatP] = boundarys(MatP, x_right, y_right)
%Apply random boundary conditions to bacteria

    jj = find(abs(MatP(:,1))>=x_right);
    MatP(jj,1) = x_right;
    MatP(jj,3) = MatP(jj,1);
    MatP(jj,2) = 2*(y_right-0.1)*rand(length(jj),1) - y_right;
    MatP(jj,4) = MatP(jj,2);
    MatP(jj,6) = 1;           %activate.
    
    jj = find(abs(MatP(:,2))>=y_right);
    MatP(jj,1) = x_right;
    MatP(jj,3) = MatP(jj,1);
    MatP(jj,2) = 2*(y_right-0.1)*rand(length(jj),1) - y_right;
    MatP(jj,4) = MatP(jj,2);
    MatP(jj,6) = 1;           %activate.
    
end
