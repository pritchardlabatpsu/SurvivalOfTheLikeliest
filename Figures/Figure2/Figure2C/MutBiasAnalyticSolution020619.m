% Analytical Solution for Two Allele Model
% Scott Leighow - 01/21/19
% 
function log2ratio = MutBiasAnalyticSolution020619(mu,M)
% 
% clear
% close all
% mu = 1e-4; % 1e-7 or 1e-6
% M = 4e6;

%% Primary Parameters

% Use CML specific parameters
% Allele A: E255V parameters
% Allele B: E255K parameters

% Fassoni LSC Dynamics Parameters
pY = 0.187/30.4;      % Fassoni 2018 - PLSC net growth (without drug) [/day]
q = -1/30.4;        % Fassoni 2018 - PLSC net growth (with drug) [/day]

turnover = 0.3;     % Komarova and Wodarz 2005 - estimated turnover
bY = pY/(1-turnover);

prob_E255V = 3e-5;      % Broad ExAC numbers
prob_E255K = 8.7e-4;

alpha_E255V = 0.0002;
alpha_E255K = 0.002;

% Birth rates [/day]
% b_S = 0.14;
% b_A = b_S;
% b_B = b_S;

b_S = bY;
b_A = b_S;
b_B = b_S;

% Net growth rates w/o drug [/day]
% r_S = 0.01;     % [/day]
% r_A = 0.0075;
% r_B = 0.005;

r_S = pY;
r_A = pY;
r_B = pY;

% Net growth rates w/ drug [/day]
% r_Sp = -0.03; 
% r_Ap = 0.0075;
% r_Bp = 0.005;
% 
% r_Ap = r_A;
% r_Bp = r_B;

r_Sp = q;
r_Ap = pY - alpha_E255V;
r_Bp = pY - alpha_E255K;

% Mutational probabilities [/mutation]
% rho_A = 0.1;
% rho_B = 0.9;

rho_A = prob_E255V/sum([prob_E255V prob_E255K]);
rho_B = prob_E255K/sum([prob_E255V prob_E255K]);

t = 0:14*365;
% t = 0:12*365;

%% Intermediate Parameters

T = 1/r_S*log(M*r_S/b_S);

alph_A = b_S^2*r_A*mu*rho_A/(r_S*b_A);
alph_B = b_S^2*r_B*mu*rho_B/(r_S*b_B);

beta_A = M*b_S*r_Ap*mu*rho_A/b_A*exp(-r_Sp*T);
beta_B = M*b_S*r_Bp*mu*rho_B/b_B*exp(-r_Sp*T);

eta_A = exp(alph_A/r_S);
eta_B = exp(alph_B/r_S);

% nu_A = eta_A*exp(-alph_A/r_S*exp(r_S*T) + beta_A/r_Sp*exp(r_Sp*T));
% nu_B = eta_B*exp(-alph_B/r_S*exp(r_S*T) + beta_B/r_Sp*exp(r_Sp*T));

nu_A = eta_A*exp(-alph_A/r_S*exp(r_S*T) + beta_A/r_Sp*exp(r_Sp*T));
nu_B = eta_B*exp(-alph_B/r_S*exp(r_S*T) + beta_B/r_Sp*exp(r_Sp*T));

%% Sensitive Population Dynamics and CDFs

% Sensitive Population
x_S(t<T) = b_S/r_S*exp(r_S*t(t<T));
x_S(t>=T) = M*exp(r_Sp*(t(t>=T)-T));

% Survival Function
S_A(t<T) =  eta_A*exp(-alph_A/r_S*exp(r_S*t(t<T)));
S_A(t>=T) = nu_A*exp(-beta_A/r_Sp*exp(r_Sp*t(t>=T)));

S_B(t<T) =  eta_B*exp(-alph_B/r_S*exp(r_S*t(t<T)));
S_B(t>=T) = nu_B*exp(-beta_B/r_Sp*exp(r_Sp*t(t>=T)));

% Cumulative Density Function
F_A = 1-S_A;
F_B = 1-S_B;

% figure
% plot(t/365,log10(x_S)/log10(M))
% hold on
% plot(t/365,[F_A; F_B]) % note: scaled so cdfs are visible
% hold off
% legend('Sensitive','F_A','F_B','Location','Northwest')

P_out1 = nu_A*nu_B;
P_out2 = (1-nu_A)*nu_B;
P_out3 = nu_A*(1-nu_B);
P_out4 = (1-nu_A)*(1-nu_B);

%% Competition

% Conditional Probability Density Function

% complete later
% f_A(t<T) = eta_A*alph_A*exp(r_S*t(t<T))*exp(-alph_A/r_S*exp(r_S*t(t<T)));
% f_A(t>=T) = nu_A*beta_A*exp(r_Sp*t(

f_A(t<T) = eta_A*alph_A*exp(r_S*t(t<T)).*exp(-alph_A/r_S*exp(r_S*t(t<T)));
f_A(t>=T) = nu_A*beta_A*exp(r_Sp*t(t>=T)).*exp(-beta_A/r_Sp*exp(r_Sp*t(t>=T)));

f_B(t<T) = eta_B*alph_B*exp(r_S*t(t<T)).*exp(-alph_B/r_S*exp(r_S*t(t<T)));
f_B(t>=T) = nu_B*beta_B*exp(r_Sp*t(t>=T)).*exp(-beta_B/r_Sp*exp(r_Sp*t(t>=T)));

fstar_A = f_A/(1-nu_A);
fstar_B = f_B/(1-nu_B);

% Times to outgrowth
Tog_A = 1/r_Ap*log(M*r_Ap/b_A);
Tog_B = 1/r_Bp*log(M*r_Bp/b_B);

delT = Tog_B - Tog_A;

% idx = 5*365:50:10*365;
% x = t(idx);
% y = x;
% [X,Y] = meshgrid(x,y);
% p = zeros(length(x));
% 
% for i = 1:(length(idx))
%     for j = 1:(length(idx))
%         p(i,j) = fstar_A(idx(i))*fstar_B(idx(j));
%     end
% end
% 
% figure
% surf(X/365,Y/365,p','FaceAlpha',0.5)
% hold on
% plot([T T]/365,[idx(1) idx(end)]/365,'--','Color',[206 28 8]/255,'LineWidth',2)
% plot([idx(1) idx(end)]/365,[T T]/365,'--','Color',[206 28 8]/255,'LineWidth',2)
% plot([idx(1) idx(end)]/365,[idx(1) idx(end)]/365-delT/365,'g--','LineWidth',2)
% % plot([t_vec(idx_vec(1)) t_vec(idx_vec(end))]/365, ([t1_vec(idx_vec(1)) t1_vec(idx_vec(end))]+delta)/365,'r','LineWidth',2)
% xlabel('t_{mut,A}')
% ylabel('t_{mut,B}')
% zlabel('prob')

% Constants
C1 = alph_A*eta_A/(1-nu_A)*alph_B*eta_B/(1-nu_B);
C2 = alph_A*eta_A/(1-nu_A)*beta_B*nu_B/(1-nu_B);
C3 = beta_A*nu_A/(1-nu_A)*alph_B*nu_B/(1-nu_B);
C4 = beta_A*nu_A/(1-nu_A)*beta_B*nu_B/(1-nu_B);

% Probability of Outcome 4.A (Integral Solution)

errbnd = 1e-4;

if delT < T     % Case 1
    
    iint_D1 = C1/(alph_B*(alph_A+alph_B*exp(-r_S*delT)))*...
        (exp(-exp(r_S*T)/r_S*(alph_A+alph_B*exp(-r_S*delT)))-...
        exp(-exp(r_S*delT)/r_S*(alph_A+alph_B*exp(-r_S*delT))))-...
        C1/(alph_A*alph_B)*exp(-alph_B/r_S)*(exp(-alph_A/r_S*exp(r_S*T))-...
        exp(-alph_A/r_S*exp(r_S*delT)));
    
%     k = 0;
%     s_k = r_Sp/r_S*(k+1);
%     v1 = alph_B/r_S*exp(r_S*(T-delT));
%     v2 = alph_B/r_S*exp(r_S*T);
%     sum_k = (-beta_A/r_Sp)^k/factorial(k)*(r_S/alph_B*exp(r_S*delT))^(k*r_Sp/r_S)*...
%             ((1-igamma(s_k,v2)/gamma(s_k))-(1-igamma(s_k,v1)/gamma(s_k)));
%     kterm = NaN;
%     while ~(abs(kterm) < errbnd)
%         k = k+1;
%         s_k = r_Sp/r_S*(k+1);
%         save1(k) = sum_k;
%         kterm = (-beta_A/r_Sp)^k/factorial(k)*(r_S/alph_B*exp(r_S*delT))^(k*r_Sp/r_S)*...
%             ((1-igamma(s_k,v2)/gamma(s_k))-(1-igamma(s_k,v1)/gamma(s_k)));
%         sum_k = sum_k + kterm;
%     end
%     
%     
    % alt
    K = -C3/(alph_B*r_S)*(r_S/alph_B*exp(r_S*delT))^(r_Sp/r_S);
    k = 0;
    v1 = alph_B/r_S*exp(r_S*(T-delT));
    v2 = alph_B/r_S*exp(r_S*T);
    dv = (v2-v1)/1000;
    v = v1:dv:v2;
    s_k = r_Sp/r_S*(k+1);
    fv_k = exp(-v).*v.^(s_k-1);
    sum_k = K*(-beta_A/r_Sp)^k/factorial(k)*...
        (r_S/alph_B*exp(r_S*delT))^(k*r_Sp/r_S)*...
        trapz(v,fv_k);
    kterm = NaN;
    save1(k+1) = sum_k;
    while ~(abs(kterm) < errbnd)
        k = k+1;
        s_k = r_Sp/r_S*(k+1);
        fv_k = exp(-v).*v.^(s_k-1);
        kterm = K*(-beta_A/r_Sp)^k/factorial(k)*...
            (r_S/alph_B*exp(r_S*delT))^(k*r_Sp/r_S)*...
            trapz(v,fv_k);
        if isnan(kterm)
            kterm = 0;
        end
        sum_k = sum_k + kterm;
        save1(k+1) = kterm;
    end
    iint_D2_1 = sum_k;
    
    iint_D2_2 = -C3/(beta_A*alph_B)*exp(-alph_B/r_S)*...
        (exp(-beta_A/r_Sp*exp(r_Sp*(T+delT)))-...
        exp(-beta_A/r_Sp*exp(r_Sp*T)));
    
    iint_D2 = iint_D2_1 + iint_D2_2;
    % Issue with conditionality: as t gets large, calculates large value
    % for integral (Inf) although the probability in that regime is
    % infintesimally small (0).  Since probability is known to converge,
    % value must approach zero (0).
    if C3 == 0
        iint_D2 = 0;
    end
    
    iint_D3 = C3/(beta_A*alph_B)*(1-exp(-beta_A/r_Sp*exp(r_Sp*(T+delT))))*...
        (exp(-alph_B/r_S*exp(r_S*T))-exp(-alph_B/r_S));
    if C3 == 0
        iint_D3 = 0;
    end
    
    iint_D4 = C4/(beta_B*(beta_A+beta_B*exp(-r_Sp*delT)))*...
        (1-exp(-exp(r_Sp*(T+delT))/r_Sp*(beta_A-beta_B*exp(-r_Sp*delT))))-...
        C4/(beta_A*beta_B)*exp(-beta_B/r_Sp*exp(r_Sp*T))*...
        (1-exp(-beta_A/r_Sp*exp(r_Sp*(T+delT))));
    if C4 == 0
        iint_D4 = 0;
    end
    
    iint_D = sum([iint_D1 iint_D2 iint_D3 iint_D4]);
        
else            % Case 2

    K = -C3/(alph_B*r_S)*(r_S/alph_B*exp(r_S*delT))^(r_Sp/r_S);
    k = 0;
    v1 = alph_B/r_S;
    v2 = alph_B/r_S*exp(r_S*T);
    dv = (v2-v1)/1000;
    v = v1:dv:v2;
    s_k = r_Sp/r_S*(k+1);
    fv_k = exp(-v).*v.^(s_k-1);
    sum_k = K*(-beta_A/r_Sp)^k/factorial(k)*...
        (r_S/alph_B*exp(r_S*delT))^(k*r_Sp/r_S)*...
        trapz(v,fv_k);
    kterm = NaN;
    while ~(abs(kterm) < errbnd)
        k = k+1;
        s_k = r_Sp/r_S*(k+1);
        fv_k = exp(-v).*v.^(s_k-1);
        kterm = K*(-beta_A/r_Sp)^k/factorial(k)*...
            (r_S/alph_B*exp(r_S*delT))^(k*r_Sp/r_S)*...
            trapz(v,fv_k);
        if isnan(kterm)
            kterm = 0;
        end
        sum_k = sum_k + kterm;
    end
    iint_D1_1 = sum_k;
    
    iint_D1_2 = C3/(beta_A*alph_B)*exp(-alph_B/r_S)*...
        (exp(-beta_A/r_Sp*exp(r_Sp*(T+delT)))-exp(-beta_A/r_Sp*exp(r_Sp*delT)));
    
    iint_D1 = iint_D1_1 + iint_D1_2;
    
    iint_D2 = C3/(beta_A*alph_B)*(1-exp(-beta_A/r_Sp*exp(r_Sp*(T+delT))))*...
        (exp(-alph_B/r_S*exp(r_S*T))-exp(-alph_B/r_S));
    if C3 == 0
        iint_D2 = 0;
    end
    
    iint_D3 = C4/(beta_B*(beta_A+beta_B*exp(-r_Sp*delT)))*...
        (1-exp(-exp(r_Sp*(T+delT))/r_Sp*(beta_A+beta_B*exp(-r_Sp*delT))))-...
        C4/(beta_A*beta_B)*exp(-beta_B/r_Sp*exp(r_Sp*T))*...
        (1-exp(-beta_A/r_Sp*exp(r_Sp*(T+delT))));
    if C4 == 0
        iint_D3 = 0;
    end
    
    iint_D = sum([iint_D1 iint_D2 iint_D3]);
    
end

%% Probability of Dominance

P_out4A = 1-iint_D;
P_out4B = iint_D;

Pdom_A = P_out2 + P_out4*P_out4A;
Pdom_B = P_out3 + P_out4*P_out4B;

log2ratio = log2(Pdom_B/Pdom_A);

% end

% % Holomorphic Lower Incomplete Gamma Function
% 
% function out = holgammainc(X,A,errbnd)
% 
% kk = 0;
% 
% sum_curr2 = X^A*exp(-X)*X^kk/gamma(A+kk+1);
% diff2 = NaN;
% 
% while ~(diff2 < errbnd)
%     kk = kk + 1;
%     sum_next2 = sum_curr2 + X^A*exp(-X)*X^kk/gamma(A+kk+1);
%     diff2 = abs(sum_next2 - sum_curr2);
%     sum_curr2 = sum_next2;
% end
% 
% out = sum_curr2;
% 
% end
