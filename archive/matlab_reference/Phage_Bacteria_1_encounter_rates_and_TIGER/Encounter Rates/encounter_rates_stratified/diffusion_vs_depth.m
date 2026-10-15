%program to compute diffusion parameters at depth
function [D0,nu,tmin] = diffusion_vs_depth(Mat, mucin, ptype)

D0 = zeros(size(Mat,1),1);                                           %array of D0's.
nu = zeros(size(Mat,1),1);                                           %array of nu's.
tmin = zeros(size(Mat,1),1);                                         %array of nu's.
%T4 hoc Parameters
D = [3.64 3.14 2.63 2.54 1.37 0.74];                      %diffusion coeff.

if ptype ==1 %T4 wild
alpha = [1.0237 0.92924 0.82187 0.81797 0.91178 0.85635];   %power law exp.
tmin0 = [-0.0003 0.0006 0.0015 0.0015 0.0008 0.0012];
else %T4 hoc
alpha = [0.99 0.99 1.0 1.0 0.99 0.89];                      %power law exp.
tmin0 = [0 0 0 0 0 0];
end

k = find(mucin >= 0 & mucin < 0.2);
    D0(k) = D(1);
    nu(k) = alpha(1);
    tmin(k) = tmin0(1);
    

k = find(mucin >= 0.2 & mucin < 0.6);
    D0(k) = D(2);
    nu(k) = alpha(2);
    tmin(k) = tmin0(2);
    
k = find(mucin >= 0.6 & mucin < 1);
    D0(k) = D(3);
    nu(k) = alpha(3);
    tmin(k) = tmin0(3);
    
k = find(mucin >= 1 & mucin < 2);
    D0(k) = D(4);
    nu(k) = alpha(4);
    tmin(k) = tmin0(4);
    
k = find(mucin >= 2 & mucin < 4);
    D0(k) = D(5);
    nu(k) = alpha(5);
    tmin(k) = tmin0(5);
    
k = find(mucin >= 4 & mucin < 5);
    D0(k) = D(6);
    nu(k) = alpha(6);
    tmin(k) = tmin0(6);
end


