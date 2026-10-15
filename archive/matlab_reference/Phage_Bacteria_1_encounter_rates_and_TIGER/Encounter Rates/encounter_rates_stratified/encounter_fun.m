function [MatP, MatB, Nenc, Nact] = encounter_fun(MatP, MatB, Rcrit, x_right)

%count number of active particles.
Nact = sum(MatP(:,5)); 

%calculate distance of phage from bacteria.
dist = sqrt((MatP(:,1)-MatB(:,1)).*(MatP(:,1)-MatB(:,1)) + (MatP(:,2)-MatB(:,2)).*(MatP(:,2)-MatB(:,2)) ...
       + (MatP(:,3)-MatB(:,3)).*(MatP(:,3)-MatB(:,3)));
   
%tally encounters with active particles (MatP(n,5)==1).
n = find(dist<=Rcrit);               %particles within attachment distance.
Nenc = sum(MatP(n,5));              %Number of active particles in contact.
MatP(n,5) = 0;                                             %deactivate = 0. 
%sim_plotter(MatP(n,:),MatB,1,0.1,10,10,0.01,[]) %plot simulation space
%Reactivate phage outside of capture distance
k = find(dist>=x_right);             %particles within attachment distance.
MatP(k,5) = 1;                                             %deactivate = 0.




