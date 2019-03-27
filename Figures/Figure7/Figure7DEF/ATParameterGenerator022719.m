%%% Resistance Profile Generator
%%% Scott Leighow - 02/24/19

clear
close all

% Randomly generates n resistance variants with allele-specific mutational
% probabilities and degrees of resistance

rng(6);

n = 10; % number of resistance variants

bias = rand(n,1);
rho = bias/sum(bias);

% Determine distribution of imatinib kill rates

IC50 = [187.5929659 143.8211099 816.799956 1829.415151 958.1929337 1055.213218 3929.957463 1081.976601 5886.708288 236.6864932 5885.450705 428.5943185 539.3844997 443.6960523 806.2694935 394.0003883 511.5376211 1482.369318 399.4440865 256.3551215];
hill = [2.605207004 1.618652833 2.338675724 1.863416542 3.00254595 1.855040394 2.03642185 2.627211788 6.699748393 1.781976566 7.49614925 2.194352549 3.445397058 1.943983667 3.703739826 2.012923777 3.017226868 2.453026827 1.264695906 2.004005035];

dose = 444;     % [nM] Rivera et al
t_assay = 3;    % [days] length of IC50
BaF3ng = 0.700; % [/day] net growth BaF3s

viab = 1./(1+(dose./IC50).^hill);
alpha_imat = -log(viab)/t_assay;

rel_alpha_imat = alpha_imat/BaF3ng;

histogram(rel_alpha_imat(2:end),5); % first IC50/hill is WT
% resistance alphas roughly exponentially distributed
% choose simulated alphas to reflect

rel_alpha = exprnd(expfit(rel_alpha_imat),n,1);

plot(rho,rel_alpha,'.')

csvwrite('ATParameters022719.csv',[rho rel_alpha])