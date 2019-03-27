%%% Loop through ATSims
%%% Scott Leighow - 02/25/19

clear

pars = csvread('ATParameters022719.csv');

log_treat = 9:12;
log_remain = 4:8;

for i = 1:length(log_treat)
    pop_treat = 10^log_treat(i);
    for j = 1:length(log_remain)
        pop_remain = 10^log_remain(j);
        
        simout = ATSims022719(pop_treat,pop_remain,pars);
    
        filename = strcat('ATSimsT',num2str(log_treat(i),'%02.f'),...
            'R',num2str(log_remain(j),'%02.f'),'_022719.csv');
        csvwrite(filename,simout);
    
    end
end