function [MatP] = boundarys_hydro(MatP, x_right, y_right, z_right, Rp, D_para_z,D_para_y,Rcrit)
%Apply random boundary conditions to bacteria

    jj = find(abs(MatP(:,1))>=x_right);
    %
    MatP(jj,1) = x_right-Rp;%x_right-Rp; %MatP(jj,1) - (2*(x_right))*MatP(jj,1)./abs(MatP(jj,1));                                      
    MatP(jj,2) = D_para_y*(2*y_right*rand(length(jj),1) - y_right); %Reflective BC's: D_para_y*(2*(y_right-Rp)*rand(length(jj),1) - y_right);  %rand y-dir.
    MatP(jj,3) = D_para_z*(2*z_right*rand(length(jj),1) - z_right); %Reflective BC's: D_para_z*(2*(z_right-Rp)*rand(length(jj),1) - z_right);  %rand z-dir.
    MatP(jj,7) = MatP(jj,1);
    MatP(jj,8) = MatP(jj,2);
    MatP(jj,9) = MatP(jj,3);
    
    %initial x-dir position 
    MatP(jj,10)=MatP(jj,1);
    %initial y-dir position 
    MatP(jj,11)=MatP(jj,2);
    %initial z-dir position 
    MatP(jj,12)=MatP(jj,3);
    
    MatP(jj,5) = 1;                                              %activate.
    %}
    %MatP(jj,:)=[];
    
    jj = find(abs(MatP(:,2))>=y_right);
    %
    MatP(jj,1) = MatP(jj,1); %MatP(jj,1) - (2*(x_right))*MatP(jj,1)./abs(MatP(jj,1));                                      
    MatP(jj,2) = D_para_y*(MatP(jj,2) - (2*(y_right))*MatP(jj,2)./abs(MatP(jj,2))); %Reflective BC's: D_para_y*(2*(y_right-Rp)*rand(length(jj),1) - y_right);  %rand y-dir.
    MatP(jj,3) = D_para_z*MatP(jj,3); %Reflective BC's: D_para_z*(2*(z_right-Rp)*rand(length(jj),1) - z_right);  %rand z-dir.
    MatP(jj,7) = MatP(jj,1);
    MatP(jj,8) = MatP(jj,2);
    MatP(jj,9) = MatP(jj,3);
    
    %initial x-dir position 
    MatP(jj,10)=MatP(jj,1);
    %initial y-dir position 
    MatP(jj,11)=MatP(jj,2);
    %initial z-dir position 
    MatP(jj,12)=MatP(jj,3);
    
    MatP(jj,5) = 1;                                              %activate.
    %}
    
    %MatP(jj,:)=[];
    jj = find(abs(MatP(:,3))>=z_right);
    %
    MatP(jj,1) = MatP(jj,1);                                       
    MatP(jj,2) = D_para_y*(MatP(jj,2)); 
    MatP(jj,3) = D_para_z*(MatP(jj,3) - (2*(z_right))*MatP(jj,3)./abs(MatP(jj,3))); 
    
    MatP(jj,7) = MatP(jj,1);
    MatP(jj,8) = MatP(jj,2);
    MatP(jj,9) = MatP(jj,3);
    
    %initial x-dir position 
    MatP(jj,10)=MatP(jj,1);
    %initial y-dir position 
    MatP(jj,11)=MatP(jj,2);
    %initial z-dir position 
    MatP(jj,12)=MatP(jj,3);
    
    MatP(jj,5) = 1;  
    %}
    
end
