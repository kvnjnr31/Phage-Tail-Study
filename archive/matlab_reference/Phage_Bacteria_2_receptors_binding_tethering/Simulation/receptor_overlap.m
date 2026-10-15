%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Program to check for receptor overlaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output] = receptor_overlap(MatB, K, Rr)
d=2*Rr;
ii = 1;

while ((d >= 2*Rr)&&(ii < K))
    
    d = sqrt((MatB(K,1)-MatB(ii,1))^2 + (MatB(K,2)-MatB(ii,2))^2 ...
        + (MatB(K,3)-MatB(ii,3))^2);
    ii = ii +1;
end
    
if d >= 2*Rr
    output = 0;    %no overlaps
else
    output = 1;    %overlaps
end

