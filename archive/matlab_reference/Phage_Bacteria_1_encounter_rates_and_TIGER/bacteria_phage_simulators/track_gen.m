function [track] = track_gen(Mat1, D0, J, nu_k, x_ic, y_ic, sim_length, dt)

%track generator

    for j =1:J
        x_old = x_ic;                                                  %x-i.c.                   
        y_old = y_ic;                                                  %y-i.c.
        Time =0;                                                 %time i.c.
        wt=0;                                            %waiting time i.c.
        
        %solve finite diff equation:
        for n=1:sim_length

            if Time >= wt

                wx = rand - 0.5;
                wy = rand - 0.5;
               
                %Langevin x-dir.
                x_new = x_old  + sqrt(24*D0*dt)*wx;
                %Langevin y-dir.
                y_new = y_old  + sqrt(24*D0*dt)*wy;   
                
                x_old=x_new;
                y_old=y_new;
                
                %generate waiting time.  
                U = rand;
                tau_c = (dt)*((gamma(1-nu_k)*U)^(-1/(nu_k)));
                %tau_c = (dt)*((U)^(-1/(nu_k)));         
                wt = Time + tau_c;                           %waiting time.
           end

        Mat1(j,(2*n-1):(2*n))=[x_old y_old];
        Time = Time + dt;
        
        end
        
    end
    
    track = Mat1;