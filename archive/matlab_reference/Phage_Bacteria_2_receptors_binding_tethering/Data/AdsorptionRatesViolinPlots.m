dbstop if error

filename = 'C:\Users\kevin\Documents\MATLAB\VirusProject\Data\PhageBacteriaFinal.xlsx';
M = xlsread(filename);

%----------------------adsorption vs. bacteria family----------------------
%Bacteria Families
Ecoli = []; Synechococcus = []; Aeruginosa = [];
for k = 1:3
    indx = find(M(:,1)==k); 
    
    if (k==1)
    %Ecoli
    Ecoli = [Ecoli; M(indx,2)]; 
    elseif (k==2)
    %Synechococcus
    Synechococcus = [Synechococcus; M(indx,2)];
    elseif (k==3)
    %Aeruginosa
    Aeruginosa = [Aeruginosa; M(indx,2)];
    end
    
end

%summary statistics of data
BacteriaAdsorbtionData.A = log10(Ecoli);
BacteriaAdsorbtionData.B = [];
BacteriaAdsorbtionData.C = log10(Aeruginosa);
BacteriaAdsorbtionData.D = [];
BacteriaAdsorbtionData.E = log10(Synechococcus);

figure(1)
%subplot(1,2,2)
violinplot(BacteriaAdsorbtionData);
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
box on
set(gca,'XTickLabel',{'E. coli','','Ps. Aeruginosa','','Synechococcus spp.'},'FontSize',18)
ax = gca;
yticks([-9 -8 -7 -6 -5])
yticklabels({'10^{-9}','10^{-8}','10^{-7}','10^{-6}','10^{-5}'})
%ax.YAxis.Scale ="log";
ax.LineWidth=2;
ylim([-9.19 -4.7])
xlim([0 6])
xtickangle(45)
%axis square
%title('Bacteria-Phage Adsorption Rates')
ylabel('Adsorption Rate (mL/hr)')

