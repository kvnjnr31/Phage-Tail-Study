function [MatP, MatB, Nenc, Nact] = encounter_fun(MatP, MatB, Rcrit)

%count number of active particles.
Nact = sum(MatP(:,6));

%find number of encounters.
dist = sqrt(((MatP(:,1)-MatB(1,1)).^2 + (MatP(:,2)-MatB(1,2)).^2));
kk = find(dist<=Rcrit);
Nenc = sum(MatP(kk,6));   %conditional encounter index:  MatP(kk,6).
MatP(kk,6) = 0;           %deactivate. 