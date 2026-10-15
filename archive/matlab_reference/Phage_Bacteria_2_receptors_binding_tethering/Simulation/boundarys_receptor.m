function [MatP] = boundarys_receptor(MatP, x_right, y_right, z_right, Rp)
%Apply random boundary conditions to bacteria

jj = find(abs(MatP(:,1))>=x_right);

if~isempty(jj)
    MatP(jj,1) = -(x_right-Rp)*MatP(jj,1)./abs(MatP(jj,1));                                      
    %MatP(jj,2) = 
    %MatP(jj,3) = 
    MatP(jj,4) = 1 + 0*MatP(jj,1);
    MatP(jj,5) = 0*MatP(jj,2);
    MatP(jj,6) = 0*MatP(jj,3);
    
    for kk = 1:length(jj)
        if MatP(kk,8) == 0
        %tail mutation
        e_tail = 2*0.01*rand - 0.01;
        mut = MatP(jj(kk),10) + e_tail;%*MatP(jj(kk),10); 
        if (mut>1)
           MatP(jj(kk),10) = 1; 
        elseif(mut<0)
           MatP(jj(kk),10) = 0;
        else
            MatP(jj(kk),10) = mut;
        end

        %head mutation
        e_head = 2*rand - 1;
        mut = MatP(jj(kk),11) + 0*e_head*MatP(jj(kk),11); 
        if (mut>1)
           MatP(jj(kk),11) = 1; 
        elseif(mut<0)
           MatP(jj(kk),11) = 0;
        else
            MatP(jj(kk),11) = mut;
        end
        %reactivate
        MatP(kk,8) = 1;
        end
    end      
end
    
    
jj = find(abs(MatP(:,2))>=y_right);

if~isempty(jj)
%MatP(jj,1) =                                       
MatP(jj,2) = -(y_right-Rp)*MatP(jj,2)./abs(MatP(jj,2));
%MatP(jj,3) = 
MatP(jj,4) = 0*MatP(jj,1);
MatP(jj,5) = 1+0*MatP(jj,2);
MatP(jj,6) = 0*MatP(jj,3);
    
    for kk = 1:length(jj)
        if MatP(kk,8) == 0
        %tail mutation
        e_tail = 2*0.01*rand - 0.01;
        mut = MatP(jj(kk),10) + e_tail; %*MatP(jj(kk),10); 
        if (mut>1)
           MatP(jj(kk),10) = 1; 
        elseif(mut<0)
           MatP(jj(kk),10) = 0;
        else
            MatP(jj(kk),10) = mut;
        end
        
        %head mutation
        e_head = 2*rand - 1;
        mut = MatP(jj(kk),11) + 0*e_head*MatP(jj(kk),11); 
        if (mut>0.5)
           MatP(jj(kk),11) = 0.5; 
        elseif(mut<0.05)
           MatP(jj(kk),11) = 0.05;
        else
            MatP(jj(kk),11) = mut;
        end
        
        %reactivate
        MatP(kk,8) = 1;
        end
    end  
end
    
jj = find(abs(MatP(:,3))>=z_right);
%
if~isempty(jj)
%MatP(jj,1) =                                        
%MatP(jj,2) = 
MatP(jj,3) = -(z_right-Rp)*MatP(jj,3)./abs(MatP(jj,3)); 

MatP(jj,4) = 0*MatP(jj,1);
MatP(jj,5) = 0*MatP(jj,2);
MatP(jj,6) = 1+0*MatP(jj,3);
    
    for kk = 1:length(jj)
        if MatP(kk,8) == 0
        %tail mutation
        e_tail = 2*0.01*rand - 0.01;
        mut = MatP(jj(kk),10) + e_tail;%*MatP(jj(kk),10); 
        if (mut>1)
           MatP(jj(kk),10) = 1; 
        elseif(mut<0)
           MatP(jj(kk),10) = 0;
        else
            MatP(jj(kk),10) = mut;
        end
        %head mutation
        e_head = 2*rand - 1;
        mut = MatP(jj(kk),11) + 0*e_head*MatP(jj(kk),11); 
        if (mut>0.5)
           MatP(jj(kk),11) = 0.5; 
        elseif(mut<0.05)
           MatP(jj(kk),11) = 0.05;
        else
            MatP(jj(kk),11) = mut;
        end
        %reactivate
        MatP(kk,8) = 1;
        end
    end  
end

  
end
