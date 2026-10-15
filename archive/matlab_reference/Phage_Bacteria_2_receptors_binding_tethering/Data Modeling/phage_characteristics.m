%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Program to Plot Distros of Phage Head and Tail Sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kBT = 300*(1.38e-23);                                  %thermal temerature.
eta = 0.001;                                           %viscocity of water.

%Myo
tails = mean([113 170 135 120 150 170 160 53.4])
heads = mean([111 85 60 60 94 64 144.5 85])

gamma0 = 6*pi*eta*(heads*(1e-9))/2;
D_myo = kBT/gamma0

figure(1)
hist(tails,3)

figure(2)
hist(heads,3)

%Sypho
tails = mean([151 180 135 980 210 100 140 355])
heads = mean([60 80 60 56 55 53 60 81])

gamma0 = 6*pi*eta*(heads*(1e-9))/2;
D_sypho = kBT/gamma0

%Podo
tails = mean([20 15 12 10 18 13 17 46])
heads = mean([60 70 63 70 59 145 60 67 60 54])

gamma0 = 6*pi*eta*(heads*(1e-9))/2;
D_Podo = kBT/gamma0
