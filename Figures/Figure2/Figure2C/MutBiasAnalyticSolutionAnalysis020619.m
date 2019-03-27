% Loop through MutBiasAnalyticSolution
% Scott Leighow - 01/22/19

clear
close all

mu_vec = logspace(-9,-6,125);
M_vec = logspace(4,12,125);

n_mu = length(mu_vec);
n_M = length(M_vec);

log2rat = NaN(n_mu,n_M);

for i = 1:n_mu
    mu = mu_vec(i);
    for j = 1:n_M
        M = M_vec(j);
        
        log2rat(i,j) = MutBiasAnalyticSolution020619(mu,M);
        
    end    
end

figure
surf(log10(mu_vec),log10(M_vec),log2rat')
xlabel('Mutation Rate')
ylabel('Population Size')

figure
heatmap = NaN(size(log2rat));
heatmap(log2rat>0)=1;
heatmap(log2rat==0)=0;
heatmap(log2rat<0)=-1;
surf(log10(mu_vec),log10(M_vec),heatmap')
xlabel('Mutation Rate')
ylabel('Population Size')
colorbar
view(0,90);

log2rat_res = reshape(log2rat,[],1);

output = [combvec(log10(mu_vec),log10(M_vec))' log2rat_res];

csvwrite('Log2Ratio.csv',output);