function [MatP, MatB, Nenc, Nact, res_times, dep_vels, N_dep] =...
         encounter_fun_hydro(MatP, MatB, res_times, Rcrit, dep_vels, Rp, Time, D_para_z,D_para_y,dt,LX)

%count number of active particles.
Nact = sum(MatP(:,5)); 

%calculate distance of phage from bacteria.
dist = sqrt(((MatP(:,1)).*(MatP(:,1)) + D_para_y*(MatP(:,2)).*(MatP(:,2)) ...
       + D_para_z*(MatP(:,3)).*(MatP(:,3))));
   
%update residence times of active particles.
n = find(dist<=Rcrit+Rp);      %find particles within residence distance.
MatP(n,6) = MatP(n,6) + dt*MatP(n,5);              %update residence times.
N_dep = 100*length(n)/10000;                           %deposited phage.

%tally encounters with active particles (MatP(n,5)==1).
n = find(dist<=Rcrit);               %particles within attachment distance.
%basin_plotter(MatP(n,11), MatP(n,12), LX)

Nenc = sum(MatP(n,5));           %conditional encounter index:  MatP(kk,5).

if(~isempty(n)&&(Time>=0))
           
j = find(MatP(n,5)==1);         %note: use actives so zeros are not stored!

%res_times = [res_times; MatP(n(j),6)];             %store residence times.             

%compute normal velocity
dist_v = [MatP(n(j),7), ... 
          D_para_y*(MatP(n(j),8)), ...
          D_para_z*(MatP(n(j),9))];
vel_v = [(MatP(n(j),1)-MatP(n(j),7))/dt, ... 
         D_para_y*(MatP(n(j),2)-MatP(n(j),8))/dt, ...
         D_para_z*(MatP(n(j),3)-MatP(n(j),9))/dt];
      
dot = vel_v(:,1).*dist_v(:,1) + D_para_y*vel_v(:,2).*dist_v(:,2) + D_para_z*vel_v(:,3).*dist_v(:,3);  %dot product at each pt.(scalar).
proj = dot./sqrt(sum((dist_v.*dist_v),2));
%dep_vels = [dep_vels; -proj];       %store deposition velocities: use -v_p.

%Event-driven detection and correction of particle positions---------------
vel_v = MatP(n,1:3)-MatP(n,7:9);                          %velocity vector.
r_old = MatP(n,7:9);                                  %old position vector.
v = sum(vel_v.*vel_v,2);                                %magnitude squared.

b = 2*sum(r_old.*vel_v,2)./v;
c = (sum(r_old.*r_old,2)-Rcrit^2)./v;
tau = -(b+sqrt(b.^2 - 4*c))/2;                          %quadratic formula.

rc = vel_v.*tau + r_old;                                   %contact points.
q = MatP(n,1:3)-rc;                                       %residual vector.
perp = -rc.*sum(rc.*q,2)./(sum(rc.*rc,2));              %correction vector.

MatP(n,1:3) = rc + q+2*perp; %updated points!
%--------------------------------------------------------------------------

end

MatP(n,5) = 0;                                             %deactivate = 0. 
MatP(n,6) = 0;                                      %reset residence times.

%reset residence times of particles outside residence distance.
n = find(dist>Rcrit+Rp);
MatP(n,6) = 0;





