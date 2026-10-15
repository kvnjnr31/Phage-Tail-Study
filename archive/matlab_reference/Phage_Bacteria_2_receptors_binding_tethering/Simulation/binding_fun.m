function [Nenc, Nact, MatP] = binding_fun(MatP, MatB, Rr)

%count number of active phage.
active = find(MatP(:,9)>0);
Nact = length(active); %number of active (unbound) phage
Nenc = 0; 

%calculate distance of phage from bacteria receptors.
for k = 1:size(MatP,1)

    if (MatP(k,8)==0) %test if unbound to any receptors
    %radius of cone of interation
    Rcone = MatP(k,10)+MatP(k,11);

    %receptor distances
    dist = sqrt((MatB(:,1)-MatP(k,1)).^2 + (MatB(:,2)-MatP(k,2)).^2 ...
    + (MatB(:,3)-MatP(k,3)).^2);        %distance of all recptors from phage k.

    indx = find(dist<Rcone);                       %receptors within tail cone.
    Pint = length(indx)*(Rr/Rcone)^2;              %probability of interaction.

        if ((rand>Pint)&&(Pint>0))                       %test for binding. 
           kk = ceil(length(indx)*rand);               %random index value.
           k_bind = indx(kk);                         %receptor that binds.
           MatP(k,8) = k_bind;                     %phage binding receptor.
           Nenc = Nenc + 1;                              %tally encounters.
        end
    end

%note: deactivation can be the prob the phage binds to three receptors
%MatP(indx,9) = 0;                                                 %detach.
end       
   
      
