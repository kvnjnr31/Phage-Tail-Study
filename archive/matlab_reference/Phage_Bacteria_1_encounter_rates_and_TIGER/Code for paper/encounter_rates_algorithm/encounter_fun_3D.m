%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Program to find encounters between bacteria and phage
%Input parameters
%MatP MatB: Phage, bacteria matrix
%Rcrit: encounter distance

%Output parameters
%MatP MatB: Phage, bacteria matrix
%Nenc: number of encounters
%Nact: number of active phage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MatP, MatB, Nenc, Nact] = encounter_fun_3D(MatP, MatB, Rcrit)

%count number of active particles.
Nact = sum(MatP(:,5));

%find number of encounters.
dist = sqrt(((MatP(:,1)-MatB(1,1)).^2 + (MatP(:,2)-MatB(1,2)).^2 + (MatP(:,3)-MatB(1,3)).^2));
kk = find(dist<=Rcrit);

Nenc = sum(MatP(kk,5));          %conditional encounter index:  MatP(kk,6).
MatP(kk,5) = 0;                                                %deactivate. 

kk = find(dist>5);
MatP(kk,5) = 1;                                                %reactivate.

