%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Algorithm to compare flux models with data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load Data
filename = 'C:\Users\kjoiner\Documents\MATLAB\adsorption_data.xlsx';
M = xlsread(filename);

%Data storage
Jball_data = [];
Jdiff_data = [];
A_data = [];

for k = 1:7
    indx = find(M(:,19)==k);
    tail_length = mean(M(indx,18));
    phage_radius = mean(M(indx,16));
    bacteria_radius = mean(M(indx,10));
    Rc = tail_length + phage_radius + bacteria_radius;
    
    bacteria_speed = mean(M(indx,11))*3600;
    bacteria_runtime = mean(M(indx,12));
    D = bacteria_speed*bacteria_speed*bacteria_runtime/3;
    
%Adsorption Rate Model
Jball = pi*Rc*Rc*bacteria_speed/(1e-6); %ballisitic
Jdiff = Rc*D; %diffusive

%Adsorption Rate Experiment
A_rate = mean(M(indx,4));

%save data
Jball_data = [Jball_data,Jball];
Jdiff_data = [Jdiff_data,Jdiff];
A_data = [A_data,A_rate];
end

figure(1)
hold on
plot(A_data./(Jball_data),'ko','MarkerSize',10)


figure(2)
plot(Jdiff_data, A_data,'ko')
%------------------------model T ~ S---------------------------------------
%linear fit
tbl = table(T',(S-mean(S))');
tbl.Properties.VariableNames = {'T', 'S'};
mdl = fitlm(tbl,'T ~ S');
beta_T=mdl.Coefficients.Estimate;  % the coefficients
