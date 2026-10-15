function [MatP] = capture_fun(MatP,Rcrit,Time,dt)

%find phage within perception distance.
dist = sqrt(((MatP(:,1)).^2 + (MatP(:,2)).^2 ...
       + (MatP(:,3)).^2));

n = find(dist>=Rcrit);
k = find(MatP(n,6)>0);

MatP(n(k),5) = Time-dt;                                  %record capture time.
MatP(n(k),6) = 0;                                              %deactivate. 

