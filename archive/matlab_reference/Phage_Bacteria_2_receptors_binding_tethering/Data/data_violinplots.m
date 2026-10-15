dbstop if error

filename = 'C:\Users\kevin\Documents\MATLAB\VirusProject\Data\bootstrap_data.xlsx';
M = xlsread(filename);

%------------------adsorption vs. phage vs. ecoli family-------------------
sipho_ecoli = []; sipho_syn = []; sipho_aeru = []; 
myo_ecoli = []; myo_syn = [];  myo_vibrio =[]; 
podo_ecoli = []; podo_aeru =[]; podo_vibrio = [];

for k = 1:9
    indx = find(M(:,1)==k); 
    
    if (k==1)||(k==2)||(k==3) %sipho (blue)
        %Ecoli
        jj = find(M(indx,2)==1);
        sipho_ecoli = [sipho_ecoli; M(indx(jj),3)]; 
        
        %Synechococcus
        jj = find(M(indx,2)==2);
        sipho_syn = [sipho_syn; M(indx(jj),3)]; 
        
        %Ps. Aeruginosa
        jj = find(M(indx,2)==3);
        sipho_aeru = [sipho_aeru; M(indx(jj),3)]; 
        
    elseif (k==4)||(k==5)||(k==6) %myo (green)
        %Ecoli
        jj = find(M(indx,2)==1);
        myo_ecoli = [myo_ecoli; M(indx(jj),3)]; 
        
        %Synechococcus
        jj = find(M(indx,2)==2);
        myo_syn = [ myo_syn; M(indx(jj),3)]; 
        
        %Vibrio
        jj = find(M(indx,2)==4);
        myo_vibrio = [myo_vibrio; M(indx(jj),3)];
        
    elseif (k==7)||(k==8)||(k==9) %podo (red)
        %Ecoli
        jj = find(M(indx,2)==1);
        podo_ecoli = [podo_ecoli; M(indx(jj),3)]; 
        
        %Ps. Aeruginosa
        jj = find(M(indx,2)==3);
        podo_aeru = [podo_aeru; M(indx(jj),3)]; 
        
        %Vibrio
        jj = find(M(indx,2)==4);
        podo_vibrio = [podo_vibrio; M(indx(jj),3)];
    end
    
end

%family adsorbtion rates
sipho_rates = [sipho_ecoli; sipho_syn; sipho_aeru];
myo_rates = [myo_ecoli; myo_syn; myo_vibrio];
podo_rates = [podo_ecoli; podo_aeru; podo_vibrio];

%----------------------------Phage Adsorption Rates------------------------
AdsorbtionData.Podoviridae = -log(podo_rates);
AdsorbtionData.Myoviridae = -log(myo_rates);
AdsorbtionData.Siphoviridae = -log(sipho_rates);

figure(2)
subplot(1,2,2)
violinplot(AdsorbtionData)
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
box on
ax = gca;
ax.YAxis.Scale ="log";
%ylim([5e-9 2e-6])
xlim([0 4])
xtickangle(45)
title('(B) Phage Adsorption Rates')
ylabel('Adsorption Rate (mL/hr)')

%-----------------------Phage-Host Adsorption Rates------------------------
PhageHostData.A = podo_ecoli;
PhageHostData.B = podo_aeru;
PhageHostData.C = podo_vibrio;
PhageHostData.D = [];
PhageHostData.E = myo_ecoli;
PhageHostData.F = myo_syn;
PhageHostData.G = myo_vibrio;
PhageHostData.H = [];
PhageHostData.I = sipho_ecoli;
PhageHostData.J = sipho_syn;
PhageHostData.K = sipho_aeru;

figure(1)
subplot(1,2,1)
violinplot(PhageHostData)
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
box on
set(gca,'XTickLabel',{'E.coli','Ps. Aeruginosa','Vibrio spp.','','E.coli','Synechococcus spp.','Vibrio spp.','','E.coli','Synechococcus spp.','Ps. Aeruginosa'},'FontSize',18)
ax = gca;
ax.YAxis.Scale ="log";
ylim([5e-9 2e-6])
xlim([0 12])
xtickangle(45)
title('(C) Phage-Host Adsorption Rates')
ylabel('Adsorption Rate (mL/hr)')

%-------------------------tail length vs. family data ---------------------
%family tail lengths
sipho_tail = [1.50e-07, 1.51e-07, 1.80e-07, 1.35e-07, 9.80e-08, 2.10e-07, 1.00e-07, 1.40e-07, 3.55e-07];
myo_tail = [1.13e-07, 1.70e-07, 1.35e-07, 1.20e-07, 1.50e-07, 1.70e-07, 1.60e-07, 5.30e-08];
podo_tail = [2.00e-08, 1.50e-08,	1.20e-08, 1.00e-08,	1.80e-08, 1.30e-08,	1.70e-08, 4.60e-08];

sipho_tail = sipho_tail'/(1e-9); myo_tail = myo_tail'/(1e-9); podo_tail = podo_tail'/(1e-9);

%family head sizes
sipho_head = 0.5*[60e-9, 60e-9, 80e-9, 60e-9 55e-9 53e-9 60e-9 81e-9];
myo_head = 0.5*[85e-9, 60e-9, 60e-9, 94e-9, 64e-9, 144.5e-9, 85.1e-9];
podo_head = 0.5*[60e-9, 70e-9, 63e-9 70e-9, 61e-9 60e-9 67e-9 60e-9];

sipho_head = sipho_head'/(1e-9); myo_head = myo_head'/(1e-9); podo_head = podo_head'/(1e-9);

TailHeadData.A = podo_tail;
TailHeadData.B = podo_head;
TailHeadData.C = [];
TailHeadData.D = myo_tail;
TailHeadData.E = myo_head;
TailHeadData.F = [];
TailHeadData.G = sipho_tail;
TailHeadData.H = sipho_head;

figure(2)
subplot(1,2,1)
violinplot(TailHeadData)
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
box on
set(gca,'XTickLabel',{'Tail','Head','','Tail','Head','','Tail','Head'},'FontSize',18)
ax = gca;
ax.YAxis.Scale ="log";
ylim([5 1e3])
xlim([0 9])
xtickangle(45)
title('(A) Phage Tail & Head Summary')
ylabel('Tail Length (nm)')

%----------------------adsorption vs. bacteria family----------------------
%family tail lengths
Ecoli = []; Synechococcus = []; Vibrio = []; Aeruginosa = [];
for k = 1:9
    indx = find(M(:,1)==k); 
    
    if (k==1)||(k==5)||(k==9)
    %Ecoli
    Ecoli = [Ecoli; M(indx,3)]; 
    elseif (k==2)||(k==4)
    %Synechococcus
    Synechococcus = [Synechococcus; M(indx,3)];
    elseif (k==6)||(k==8)
    %Vibrio
    Vibrio = [Vibrio; M(indx,3)];
    elseif (k==3)||(k==7)
    %Aeruginosa
    Aeruginosa = [Aeruginosa; M(indx,3)];
    end
    
end

%summary statistics of data
BacteriaAdsorbtionData.A = Synechococcus;
BacteriaAdsorbtionData.B = [];
BacteriaAdsorbtionData.C = Ecoli;
BacteriaAdsorbtionData.D = [];
BacteriaAdsorbtionData.E = Vibrio;
BacteriaAdsorbtionData.F = [];
BacteriaAdsorbtionData.G = Aeruginosa;

figure(1)
subplot(1,2,2)
violinplot(BacteriaAdsorbtionData)
set(gca,'FontSize',18,'LineWidth',2,'FontWeight','bold')
box on
set(gca,'XTickLabel',{'Synechococcus spp.','','E. coli','','Vibrio spp.','','Ps. Aeruginosa'},'FontSize',18)
ax = gca;
ax.YAxis.Scale ="log";
ax.LineWidth=2;
ylim([5e-9 2e-6])
xlim([-1 8])
xtickangle(45)
title('(D) Host-Phage Adsorption Rates')
ylabel('Adsorption Rate (mL/hr)')

%{
%---------------statistical modeling of predicted encounter rates----------
%Ecoli
r_ecoli = 0.5*(1.19e-6); v_ecoli = 21e-6;
%RC = r_ecoli+[sipho_tail*(1e-9);myo_tail*(1e-9);podo_tail*(1e-9)];

taildata0 = [sipho_tail*(1e-9);myo_tail*(1e-9);podo_tail*(1e-9)];
headdata0 = [sipho_head*(1e-9);myo_head*(1e-9);podo_head*(1e-9)];

for kk = 1:10000
   indx = randsample([1 2 3],1);
   if indx == 1
      headdata = sipho_head*(1e-9);
      taildata = sipho_tail*(1e-9);
   elseif indx ==2
      headdata = myo_head*(1e-9);
      taildata = myo_tail*(1e-9);
   elseif indx == 3
      headdata = podo_head*(1e-9);
      taildata = podo_tail*(1e-9);
   end
   
   head = randsample(headdata,1);
   tail = randsample(taildata,1);
   RC(kk) = r_ecoli + head + tail;
end
ecoli_rate = (1e6)*3600*pi*v_ecoli*RC.^2;
k_ecoli = median(ecoli_rate);

t_ecoli = median(RC)/v_ecoli

%Synechococcus
r_syn = 0.5*1.04e-6; v_syn = 15e-6;
RC = r_syn+[sipho_tail*(1e-9);myo_tail*(1e-9)];
taildata0 = [sipho_tail*(1e-9);myo_tail*(1e-9)];
headdata0 = [sipho_head*(1e-9);myo_head*(1e-9)];

for kk = 1:10000
   indx = randsample([1 2],1);
   if indx == 1
      headdata = sipho_head*(1e-9);
      taildata = sipho_tail*(1e-9);
   elseif indx ==2
      headdata = myo_head*(1e-9);
      taildata = myo_tail*(1e-9);
   elseif indx == 3
      headdata = podo_head*(1e-9);
      taildata = podo_tail*(1e-9);
   end
   
   head = randsample(headdata,1);
   tail = randsample(taildata,1);
   RC(kk) = r_syn + head + tail;
end
syn_rate = (1e6)*3600*pi*v_syn*RC.^2;
k_syn = median(syn_rate);
t_syn = median(RC)/v_syn

%Aeruginosa
r_areu = 0.5*1.19e-6; v_areu = 51.3e-6;
RC = r_areu+[sipho_tail*(1e-9);podo_tail*(1e-9)];
taildata0 = [sipho_tail*(1e-9);podo_tail*(1e-9)];
headdata0 = [sipho_head*(1e-9);podo_head*(1e-9)];
for kk = 1:1000000
   indx = randsample([1 3],1);
   if indx == 1
      headdata = sipho_head*(1e-9);
      taildata = sipho_tail*(1e-9);
   elseif indx ==2
      headdata = myo_head*(1e-9);
      taildata = myo_tail*(1e-9);
   elseif indx == 3
      headdata = podo_head*(1e-9);
      taildata = podo_tail*(1e-9);
   end
   head = randsample(headdata,1);
   tail = randsample(taildata,1);
   RC(kk) = r_areu + head + tail;
end
areu_rate = (1e6)*3600*pi*v_areu*RC.^2;
k_areu = median(areu_rate);
t_areu = median(RC)/v_areu

%Vibrio
r_vibrio = 0.5*0.72e-6; v_vibrio = 73.7e-6;
RC = r_vibrio+[podo_tail*(1e-9);myo_tail*(1e-9)];
taildata0 = [podo_tail*(1e-9);myo_tail*(1e-9)];
headdata0 = [podo_head*(1e-9);myo_head*(1e-9)];
for kk = 1:10000
   indx = randsample([2 3],1);
   if indx == 1
      headdata = sipho_head*(1e-9);
      taildata = sipho_tail*(1e-9);
   elseif indx ==2
      headdata = myo_head*(1e-9);
      taildata = myo_tail*(1e-9);
   elseif indx == 3
      headdata = podo_head*(1e-9);
      taildata = podo_tail*(1e-9);
   end
   head = randsample(headdata,1);
   tail = randsample(taildata,1);
   RC(kk) = r_vibrio + head + tail;
end
vibrio_rate = (1e6)*3600*pi*v_vibrio*RC.^2;
k_vibrio = median(vibrio_rate);
t_vibrio = median(RC)/v_vibrio

%------------------------test of host adsorption data----------------------
%Ecoli vs Synechococcus
[p, h, stats] = ranksum(Synechococcus,Ecoli);
if h ~=0
    disp('Ecoli vs Synechococcus')
    p
end
%Ecoli vs Aeruginosa
[p, h, stats] = ranksum(Aeruginosa,Ecoli);
if h ~=0
    disp('Ecoli vs Aeruginosa')
    p
end
%Ecoli vs Vibrio
[p, h stats] = ranksum(Vibrio,Ecoli);
if h ~=0
   disp('Ecoli vs Vibrio')
    p
end
%Synechococcus vs Aeruginosa
[p, h stats] = ranksum(Aeruginosa,Synechococcus);
if h ~=0
    disp('Synechococcus vs Aeruginosa')
    p
end
%Synechococcus vs Vibrio
[p, h stats] = ranksum(Vibrio,Synechococcus);
if h ~=0
    disp('Synechococcus vs Vibrio')
    p
end
%Aeruginosa vs Vibrio
[p, h stats] = ranksum(Vibrio,Aeruginosa);
if h ~=0
    disp('Aeruginosa vs Vibrio')
    p
end

%--------------------------------predicted rates---------------------------
rates = [median(syn_rate);   median(ecoli_rate); median(vibrio_rate); median(areu_rate)];
max_rates = [prctile(syn_rate,75);  prctile(ecoli_rate,75); prctile(vibrio_rate,75); prctile(areu_rate,75)];
min_rates = [prctile(syn_rate,25);  prctile(ecoli_rate,25); prctile(vibrio_rate,25); prctile(areu_rate,25)];
positions = [0 2 4 6];

%Others
indx = [1 2 3 4];
h2=errorbar(positions(indx),rates(indx),rates(indx)-min_rates(indx),max_rates(indx)-rates(indx),'ko','LineStyle','none','LineWidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','k','MarkerFaceColor','w');
h2.CapSize = 12;
%}



