%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Adsorption Model Bootstrap Algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'C:\Users\kjoiner\Documents\MATLAB\Virus Project\Data\bootstrap_data.xlsx';
M = xlsread(filename);

for trails = 1:2000
%----------------------------Extract Data----------------------------------
data = [];
for k = 1:10
    indx = find(M(:,1)==k);
    
    Adsorb_rate=M(indx,2);
    
    %phage tail parameter
    tail_length = [];
    for j = 1: length(indx)
    tails = M(indx(j),5:13);
    tail_indx = tails>0;
    tail_length(j,1) = datasample(tails(tail_indx),1);%M(indx,1);
    end
    
    %phage head parameter
    head_size = [];
    for j = 1: length(indx)
    heads = M(indx(j),14:23);
    head_indx = tails>0;
    head_size(j,1) = datasample(heads(head_indx),1);%M(indx,1);
    end
    
    %bacteria parameter
    dmin = M(indx,3); dmax = M(indx,4);
    bacteria_diameter = (dmax + dmin)/2;%.*rand(length(indx),1) + dmin;
    
       
    %store data
    data0 = [mean(bacteria_diameter), mean(tail_length), mean(head_size), mean(Adsorb_rate)];  
    data = [data;data0];
end

%----------------------Adsorption Rate linear fit--------------------------
B = data(:,1); T = data(:,2); H = data(:,3); A = data(:,4);
B = (B-mean(B))/std(B); T = (T-mean(T))/std(T);
H = (H-mean(H))/std(H); A = (A-mean(A))/std(A);
%{
figure(199)
plot3(B, T, A,'ko')
zlabel('Adsorption Rate (mL/hr)')
xlabel('bacteria radius (nm)')
ylabel('tail length (nm)')
%}

tbl = table(B, T, H, A);
tbl.Properties.VariableNames = {'B', 'T', 'H', 'A'};
mdl = fitlm(tbl,'A ~ B + T + H');

betas = mdl.Coefficients.Estimate;  %model coefficents
pvalues = mdl.Coefficients.pValue;  %model pvalues

R2(trails) = mdl.Rsquared.Ordinary;
tail_pvalue(trails) = pvalues(3);
bacteria_pvalue(trails) = pvalues(2);

bacteria_coeff(trails) = betas(2);
tail_coeff(trails) = betas(3);
head_coeff(trails) = betas(4);
end

figure(99)
hist(R2,20);

[counts value] = hist(R2,20);
normal_counts = counts/sum(counts); %normalization

cdf = cumsum(normal_counts);

data_mat = [R2', tail_pvalue', bacteria_pvalue'];
data_mat = sortrows(data_mat);

figure(1)
hold on
yyaxis right
plot(data_mat(:,1), data_mat(:,3),'r', data_mat(:,1), data_mat(:,2),'m','LineWidth',2)
yyaxis left
plot(value, 1-cdf,'k','LineWidth',2)

figure(2)
hist(bacteria_coeff,20)

figure(3)
hist(head_coeff,20)

figure(4)
hist(tail_coeff,20)

aa=1;
%{
[Bx,Tx] = meshgrid(linspace(min(0.8*B),max(1.2*B),100),linspace(min(0.8*T),max(1.2*T),100));
fit = betas(1) + betas(2)*Bx + betas(3)*Tx;

figure(2)
hold on
surf(Bx, Tx, fit,'EdgeColor','none','FaceAlpha',0.75)
colormap('gray');
%view(rot_indx,37)

%---------------------Adsorption Rate Experimental Plots------------------- 
%sipho (blue)
plot3(data(1:3,1)/(1e-9),data(1:3,2)/(1e-9), (data(1:3,3)),'ko','Markersize',12,'Linewidth',2,'MarkerFaceColor',[0 0 1])
%myo (green)
plot3(data(4:7,1)/(1e-9),data(4:7,2)/(1e-9),(data(4:7,3)),'ko','Markersize',12,'Linewidth',2,'MarkerFaceColor',[0 1 0])
%podo (red)
plot3(data(8:10,1)/(1e-9),data(8:10,2)/(1e-9),(data(8:10,3)),'ko','Markersize',12,'Linewidth',2,'MarkerFaceColor',[1 0 0])

title('Adsorption Rate v.s. Tail Length (T=37 deg C)')
zlabel('Adsorption Rate (mL/hr)')
xlabel('bacteria radius (nm)')
ylabel('tail length (nm)')
%axis([0 0.17 0 3])
set(gca,'FontSize',22)
%}
