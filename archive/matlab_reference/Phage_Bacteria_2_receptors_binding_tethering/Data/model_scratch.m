%model scratch

filename = 'C:\Users\kjoiner\Documents\MATLAB\Virus Project\Data\adsorption_data.xlsx';
M = xlsread(filename);

%-----------------------adsorption vs. tail length-------------------------
data=[]; Jball_data = []; Jdiff_data = []; sipho = []; myo = []; podo = [];

for k = 1:10
    indx = find(M(:,19)==k);

%--------------------------Adsorption Rate Models--------------------------
   %Parameters
    tail_length = mean(M(indx,18)); 
    phage_radius = mean(M(indx,16)); 
    bacteria_radius = mean(M(indx,10));
    Rc = tail_length + phage_radius + bacteria_radius;
    bacteria_speed = mean(M(indx,11));
    bacteria_runtime = mean(M(indx,12));
    D = bacteria_speed*bacteria_speed*bacteria_runtime/3;
    
    %Models
    Jball = (36e8)*pi*Rc*Rc*bacteria_speed; %ballisitic
    
    if isnan(D)
        Jdiff = Jball;
    else
    Jdiff = (36e8)*Rc*D; %diffusive
    end
    
    %Save data
    data0 = [Jball];
    Jball_data = [Jball_data; data0];

    data0 = [Jdiff];
    Jdiff_data = [Jdiff_data; data0];

end