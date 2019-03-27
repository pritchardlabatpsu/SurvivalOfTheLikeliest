%%%% Loop through LSCMutBiasGeneral

clear
close all

%% Import values

probs = [0.038343558 0.012269939 0.133435583 0.023006135 0.004601227 0.055214724 0.133435583 0.004601227 0.038343558 0.182515337 0.076687117 0.055214724 0.038343558 0.009202454 0.003067485 0.009202454 0.010736196 0.038343558 0.133435583];

%% Loop through simulations

drug_vec = {'Imatinib' 'Nilotinib' 'Dasatinib' 'Bosutinib'};
mutbias_vec = [0 1];

% drug_vec = {'Imatinib'};
% mutbias_vec = [0 1];
% pharma_vec = 1;

for i = 1:length(drug_vec)

    drug = string(drug_vec(i));
    
    alpha_drug = csvread(strcat(drug,'AlphaGeneratorResults011619.csv'));
    
    for j = 1:length(mutbias_vec)
       
        mutbias = mutbias_vec(j);
        
        if mutbias
            biasR = probs/sum(probs);
            biasstr = 'WBias';
        else
            biasR = ones(1,length(probs))/length(probs);
            biasstr = 'WoBias';
        end

        output = LSCMutBiasGeneral_SL011619(alpha_drug,biasR);

        filename = strcat('LSCMutBias',drug,biasstr,'011619.csv');
        csvwrite(filename,output);
            
        
    end 
    
end
