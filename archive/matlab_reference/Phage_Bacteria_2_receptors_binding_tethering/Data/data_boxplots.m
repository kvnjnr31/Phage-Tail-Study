dbstop if error

filename = 'C:\Users\kjoiner\Documents\MATLAB\Virus Project\Data\bootstrap_data.xlsx';
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

%statistical test of family adsorbtion rates
%{
disp('sipho vs myo')
[p, h, stats] = ranksum(sipho_rates,myo_rates)
disp('sipho vs podo')
[p, h, stats] = ranksum(sipho_rates,podo_rates)
disp('myo vs podo')
[p, h, stats] = ranksum(myo_rates,podo_rates)
%}
%-----------------------summary statistics of data-------------------------
figure(1)
title('Adsorption Rate Summary Statistics')
subplot(1,2,1)
colors = [1 0 0 ; 0 1 0; 0 0 1]; %sipho (blue), myo (green), podo (red)

%tag groups
group = [ones(size(podo_ecoli)); 2*ones(size(podo_aeru)); 3*ones(size(podo_vibrio));
    4*ones(size(myo_ecoli)); 5*ones(size(myo_syn)); 6*ones(size(myo_vibrio));
    7*ones(size(sipho_ecoli)); 8*ones(size(sipho_syn)); 9*ones(size(sipho_aeru))];
positions = [1 1.25 1.50 2 2.25 2.50 3 3.25 3.50];

h = boxplot([podo_ecoli; podo_aeru; podo_vibrio; myo_ecoli; myo_syn; myo_vibrio; sipho_ecoli; sipho_syn; sipho_aeru], group,'positions', positions);

set(h,{'linew'},{2},'Color', 'k')
hh = findobj(gca,'Tag','Box');

for j=1:length(positions)
    patch(get(hh(j),'XData'),get(hh(j),'YData'),colors(ceil(j/3),:),'FaceAlpha',.5);
end

set(gca,'XTickLabel',{'E.coli','Ps. Aeruginosa','Vibrio spp.','E.coli','Synechococcus spp.','Vibrio spp.','E.coli','Synechococcus spp.','Ps. Aeruginosa'},'FontSize',18)
set(h,{'linew'},{2})
ax = gca;
ax.YAxis.Scale ="log";
ax.LineWidth=2;
ylim([5e-9 2e-6])
xlim([0.75 3.75])

c = get(gca, 'Children');
hleg1 = legend([c(1),c(4),c(7)], 'Podoviridae', 'Myoviridae', 'Siphoviridae' );

xtickangle(45)
title('(C) Phage-Host Adsorption Rates')
ylabel('Adsorption Rate (mL/hr)')

%----------------------------End Figure(1)---------------------------------

%stats of adsorption data
%{
podo_ecoli_med =median([podo_ecoli])
podo_aeru_med = median([podo_aeru])
podo_vibrio_med = median([podo_vibrio])
myo_ecoli_med = median([myo_ecoli])
myo_syn_med = median([myo_syn])
myo_vibrio_med = median([myo_vibrio])
sipho_ecoli_med = median([sipho_ecoli])
sipho_syn_med = median([sipho_syn])
sipho_aeru_med = median([sipho_aeru])
%}

%stats test of E.coli adsorption data
disp('sipho vs myo')
[p, h, stats] = ranksum(sipho_ecoli,myo_ecoli);
disp('sipho vs podo')
[p, h, stats] = ranksum(sipho_ecoli,podo_ecoli);
disp('myo vs podo')
[p, h, stats] = ranksum(myo_ecoli,podo_ecoli);
%-------------------------tail length vs. family---------------------------
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

%summary statistics of head tail data
figure(2)
subplot(1,2,1)
colors = [1 0 0 ;0 1 0; 0 0 1]; %sipho (blue), myo (green), podo (red)

group = [ones(size(podo_tail)); 2*ones(size(podo_head));
    3*ones(size(myo_tail)); 4*ones(size(myo_head));
    5*ones(size(sipho_tail)); 6*ones(size(sipho_head))];
positions = [1 1.25 1.62 1.87 2.24 2.49];

h = boxplot([podo_tail; podo_head; myo_tail; myo_head; sipho_tail; sipho_head;], group,'positions', positions);

set(h,{'linew'},{2},'Color', 'k')
hh = findobj(gca,'Tag','Box');

for j=1:length(positions)
    patch(get(hh(j),'XData'),get(hh(j),'YData'),colors(ceil(j/2),:),'FaceAlpha',.5);
end

set(gca,'XTickLabel',{'Tail','Head','Tail','Head','Tail','Head'},'FontSize',18)
ax = gca;
ax.YAxis.Scale ="log";
ax.LineWidth=2;
xtickangle(45)
ylim([5 1e3])
title('(A) Phage Tail & Head Summary')
ylabel('Tail Length (nm)')

subplot(1,2,2)
colors = [1 0 0 ;0 1 0; 0 0 1]; %sipho (blue), myo (green), podo (red)
group = [ ones(size(podo_rates)); 2 * ones(size(myo_rates)); 3 * ones(size(sipho_rates))];
h = boxplot([podo_rates; myo_rates; sipho_rates],group);

set(h,{'linew'},{2},'Color', 'k')
hh = findobj(gca,'Tag','Box');

for j=1:3
    patch(get(hh(j),'XData'),get(hh(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

c = get(gca, 'Children');
hleg1 = legend([c(1),c(2),c(3)], 'Podoviridae', 'Myoviridae', 'Siphoviridae' );

set(gca,'XTickLabel',{'Podoviridae','Myoviridae','Siphoviridae'},'FontSize',18)
ax = gca;
ax.YAxis.Scale ="log";
ax.LineWidth=2;
xtickangle(45)
ylim([5e-9 2e-6])
title('(B) Phage Adsorption Rates')
ylabel('Adsorption Rate (mL/hr)')
%----------------------------End Figure(2)---------------------------------

%statistical test of family tail lengths
disp('sipho vs myo')
[p, h, stats] = ranksum(sipho_tail,myo_tail);

disp('sipho vs podo')
[p, h stats] = ranksum(sipho_tail,podo_tail);

disp('myo vs podo')
[p, h, stats] = ranksum(myo_tail,podo_tail);


tail_median_sipho = median(sipho_tail);
tail_median_myo = median(myo_tail);
tail_median_podo = median(podo_tail);

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
figure(1)
subplot(1,2,2)
hold on
colors = [0.5 0.5 0.5 ;0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5]; 
group = [ ones(size(Synechococcus)); 2 * ones(size(Ecoli)); 3 * ones(size(Vibrio)); 4 *ones(size(Aeruginosa)) ];
A = 1; B = 3.75; C = 6.5; D = 9.25;
positions = [A B C D];
h = boxplot([Synechococcus;  Ecoli; Vibrio; Aeruginosa],group,'positions',positions);

set(h,{'linew'},{2},'Color', 'k')
hh = findobj(gca,'Tag','Box');

for j=1:4
    patch(get(hh(j),'XData'),get(hh(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

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
positions = [A B C D]-1.25;

%Others
indx = [1 2 3 4];
h2=errorbar(positions(indx),rates(indx),rates(indx)-min_rates(indx),max_rates(indx)-rates(indx),'ko','LineStyle','none','LineWidth',2,'MarkerSize',10,...
    'MarkerEdgeColor','k','MarkerFaceColor','w');
h2.CapSize = 12;

set(gca,'XTickLabel',{'Synechococcus spp.','E. coli','Vibrio spp.','Ps. Aeruginosa'},'FontSize',18)
ax = gca;
ax.YAxis.Scale ="log";
ax.LineWidth=2;
ylim([5e-9 2e-6])
xlim([-1 10.25])
xtickangle(45)
title('(D) Host-Phage Adsorption Rates')
%ylabel('Adsorption Rate (mL/hr)')




