
% Imatinib Pharmacokinetic Model Parameterization

clear
close all

rng(50);

%% Parameters
tic;
n = 15; % number of samples drawn for each generated parameter

%% Drug Data

drug = 'Imatinib'; % ['Imatinib' 'Dasatinib' 'Nilotinib']

drug_vec = {'Imatinib' 'Dasatinib' 'Nilotinib'};

for d = 1:length(drug_vec)
    
    rng(50);
    
    drug = drug_vec(d);
    
    if strcmp(drug,'Imatinib')

        % Raw values from Peng 2005

        t_up = 3/24;
        t_downhalf = 18/24;
        t_dose = 24/24;

        MW = 493.603;     % molecular weight [g/mol]
        SF = 3377/444;    % Rivera poster shift factor

        CI = 0.95;

        Cmax = 2.6/(MW*SF)*1e6;     % [nM]
        Cmax_sd = 0.8/(MW*SF)*1e6/norminv(0.5+CI/2); 

        Cmin = 1.2/(MW*SF)*1e6; 
        Cmin_sd = 0.8/(MW*SF)*1e6/norminv(0.5+CI/2);

        Cmax_vec = normrnd(Cmax,Cmax_sd,n,1);
        Cmin_vec = normrnd(Cmin,Cmin_sd,n,1);
        
        IC50 = [187.5929659 143.8211099 816.799956 1829.415151 958.1929337 1055.213218 3929.957463 1081.976601 5886.708288 236.6864932 5885.450705 428.5943185 539.3844997 443.6960523 806.2694935 394.0003883 511.5376211 1482.369318 399.4440865 256.3551215];
        hill = [2.605207004 1.618652833 2.338675724 1.863416542 3.00254595 1.855040394 2.03642185 2.627211788 6.699748393 1.781976566 7.49614925 2.194352549 3.445397058 1.943983667 3.703739826 2.012923777 3.017226868 2.453026827 1.264695906 2.004005035];

    elseif strcmp(drug,'Dasatinib')

        SF = 27/11; % Rivera

        % FDA access data - Sprycel NDA 22-072
        t_up = 3/24;
        t_downhalf = 4/24;
        t_dose = 24/24;

        % Raw values from Rousselot 2010

        CI = 0.95;

        % Assume lognormal distribution
        Cmin_med = 2.1/SF;
        Cmin_lo = 0.2/SF;
        Cmin_hi = 18.7/SF;

        Cmin_sd1 = (log10(Cmin_lo)-log10(Cmin_med))/norminv(.5-CI/2);
        Cmin_sd2 = (log10(Cmin_hi)-log10(Cmin_med))/norminv(.5+CI/2);
        Cmin_sd = mean([Cmin_sd1 Cmin_sd2]);

        Cmax_med = 107.6/SF;
        Cmax_lo = 20.5/SF;
        Cmax_hi = 353/SF;

        Cmax_sd1 = (log10(Cmax_lo)-log10(Cmax_med))/norminv(.5-CI/2);
        Cmax_sd2 = (log10(Cmax_hi)-log10(Cmax_med))/norminv(.5+CI/2);
        Cmax_sd = mean([Cmax_sd1 Cmax_sd2]);

        Cmax_vec = 10.^(normrnd(log10(Cmax_med),Cmax_sd,n,1));
        Cmin_vec = 10.^(normrnd(log10(Cmin_med),Cmin_sd,n,1));
        
        IC50 = [2.411092212 3.355991196 6.95860279 8.072449156 3.500165403 1.35140484 5.875668996 5.531157242 8.609195389 1.624658866 9550.906064 13.88002198 3.880493595 2.686224945 3.270415087 3.645290721 3.708111541 2.711545335 12.80454711 2.239436811];
        hill = [2.843638349 1.996274662 1.163728817 1.663090374 1.605442651 2.042049628 2.443078153 2.332149484 1.804567259 2.118828791 3.976962956 2.289360733 3.698605597 1.997628663 2.744691504 2.148451942 1.926989673 1.554698092 0.62421435 1.827723367];
        
    elseif strcmp(drug,'Nilotinib')

        SF = 2754/131;
        MW = 529.5;

        % Raw values from Trent 2011

        t_up = 3/24;
        t_downhalf = 17/24;
        t_dose = 12/24;

        Cmax = 1.644/(MW*SF)*1e6;
        Cmax_sd = .828/(MW*SF)*1e6;

        % estimate Cmin
        b_est = log(2)/t_downhalf;
        Cmin = Cmax*exp(-b_est*t_dose);
        Cmin_sd = Cmax_sd;

        Cmax_vec = normrnd(Cmax,Cmax_sd,n,1);
        Cmin_vec = normrnd(Cmin,Cmin_sd,n,1);
        
        IC50 = [15.7250575 26.33205467 54.83601896 121.5617847 44.14997012 50.44671031 720.0781776 92.42916917 576.8920447 14.95905801 45177.01919 48.21748905 20.22680974 25.90526803 122.0577871 103.8235483 107.9199748 122.4738536 70.38838611 22.62561518];
        hill = [2.117926267 2.426923242 1.309207941 1.106411692 1.629186724 1.856100477 2.470549328 2.125789995 2.480999438 1.437378229 3.409430106 1.702166415 1.855711211 1.828055724 2.435492178 1.540404962 2.230879252 1.842694962 0.734725263 1.84623322];
        
    end

    % Distribution of WBC count at detection (numbers from Qin Medicine 2016;
    % lognormal dist according to Stein 2011)
    WBC_med = 98e9;   % /L    Qin Medicine 2016: Combination of WBC Countat Presentation...
    WBC_lo = 6.8e9;
    WBC_hi = 683.5e9;

    vol_blood = 10;     % L

    % Estimate standard deviation (assume that range of WBC counts observed in 
    % Qin (N=362) fall within 95% CI)
    WBC_CI = .95;    
    WBC_sd1 = (log10(WBC_lo)-log10(WBC_med))/norminv(.5-WBC_CI/2);
    WBC_sd2 = (log10(WBC_hi)-log10(WBC_med))/norminv(.5+WBC_CI/2);

    WBC_sd = mean([WBC_sd1 WBC_sd2]);

    % Distribution of detection sizes (follows lognormal dist according to
    % Stein et al.)
    log_WBCcount = normrnd(log10(WBC_med),WBC_sd,n,1);
    WBC_count = 10.^log_WBCcount;

    %% Simulate Patients

    options = optimset('Display','off');

    tend = 365*2/t_dose;
    tstep = 0.01;

    nvars = length(IC50);

    t_assay = 3;

    BaF3_ng = 0.700; % our BaF3 count data [/day]

    % Adherence distribution - Marin 2010

    adher_range = [1 .975 .925 .875 .825 .8];   % degrees of adherence
    adher_npats = [36 21 7 5 6 12];             % number of patients observed in each range

    adher_range_int = 0.8:0.01:1;
    adher_npats_int = interp1(adher_range,adher_npats,adher_range_int,'spline');

    % Generate patient specific parameters

    adher_vec = randsample(adher_range_int,n,true,adher_npats_int)';
    poptreat_vec = WBC_count*vol_blood;    % population size at the start of treatment [cells]

    % Generate dosing schedule
    dose_mat = NaN(n,tend);
    for i = 1:n
        adher_i = adher_vec(i);
        dose_mat(i,:) = randsample([1 0],tend,true,[adher_i (1-adher_i)]); % does patient take imatinib today?
    end

    alpha_mat = NaN(n,n,n,nvars);

    parfor i = 1:n
        Cmax_i = Cmax_vec(i);
        alpha_mat_i = NaN(n,n,nvars);
        for j = 1:n
            Cmin_i = Cmin_vec(j);
            for k = 1:n
                adher_i = adher_vec(k);

                y0 = Cmin_i;

                % Estimate b from half-life: 0.5*y = y*exp(-b*t_downhalf)

                b_0 = -log(0.5)/t_downhalf;

                % Estimate a_0 from t_up: Cmax = Cmin*exp(a*t_up)

                a_0 = log(Cmax_i/Cmin_i)/t_up;

                % Estimate x0 from t_dose: y(t_dose)=Cmin

                x0_0 = (Cmin_i-Cmin_i*exp(-b_0*t_dose))*(a_0-b_0)/(a_0*exp(-b_0*t_dose));

                % Solve for x0, a, and b

                % conditions: 
                % y(tmax) = Cmax
                % dydt(tmax) = 0
                % y(t_dose) = Cmin

                Fx = @(x) [(x(2)*x(1)/(x(2)-x(3))+y0)*exp(-x(3)*t_up) - x(2)*x(1)/(x(2)-x(3))*exp(-x(2)*t_up) - Cmax_i
                           -x(3)*(x(2)*x(1)/(x(2)-x(3))+y0)*exp(-x(3)*t_up) + x(2)^2*x(1)/(x(2)-x(3))*exp(-x(2)*t_up)
                           (x(2)*x(1)/(x(2)-x(3))+y0)*exp(-x(3)*t_dose) - x(2)*x(1)/(x(2)-x(3))*exp(-x(2)*t_dose) - Cmin_i];

                x_init = [x0_0 a_0 b_0];

                [x,~,exitflag,~] = fsolve(Fx,x_init,options);

                % Continue if solution not found or if solution is complex (occurs for
                % extreme concentration values)
                if exitflag<=0 || any(imag(x))
                    continue
                end

                x0 = x(1);
                a = x(2);
                b = x(3);

                y = y0;
                t_j = (0:tstep:t_dose)';

                for m = 1:tend

                    y0_m = y(end);

                    dosed = dose_mat(k,m);

                    if dosed
                        y_m = (a*x0/(a-b)+y0_m)*exp(-b*t_j)-a*x0/(a-b)*exp(-a*t_j);
                        y = [y;y_m(2:end)];
                    else
                        y_m = y0_m*exp(-b*t_j);
                        y = [y;y_m(2:end)];
                    end

                end

                %% Transform to Alpha_Ave Values

                alpha_ave = NaN(1,nvars);

                for m = 1:nvars
                    IC50_m = IC50(m);
                    hill_m = hill(m);
                    viab_t = 1./(1+(y/IC50_m).^hill_m);
                    alpha_t = -log(viab_t)/t_assay;
                    alpha_ave(m) = mean(alpha_t);
                end

                alpha_mat_i(j,k,:) = alpha_ave;


            end
        end    
        alpha_mat(i,:,:,:) = alpha_mat_i;
        disp(i)
    end

    %% Transform to LSC Values

    alpha_rmat = reshape(alpha_mat,[],20);

    par_mat = combvec(Cmax_vec',Cmin_vec',adher_vec',log10(poptreat_vec)')';
    comp_mat = [par_mat repmat(alpha_rmat,n,1)];

    comp_mat = comp_mat(~isnan(comp_mat(:,5)),:);

    pY = 0.18701/30.4;  % PLSC net growth (without drug) [/day]
    q = 1.00394/30.4;   % PLSC net drug kill term [/day]

    % Transform WT alphas for expected net drug kill rate
    % Scale other drugs by imatinib kill rate

    WTalphaBaF3 = comp_mat(:,5);
    resalphaBaF3 = comp_mat(:,6:end);

    WTqBaF3 = BaF3_ng - WTalphaBaF3;
    
    if strcmp(drug,'Imatinib')
        WTqBaF3imat = mean(WTqBaF3);
    end

    WTqLSC = WTqBaF3/WTqBaF3imat*(-q);
    WTalphaLSC = pY - WTqLSC;

    % Transform res alphas by growth rate of each system

    resalphaLSC = resalphaBaF3/BaF3_ng*pY;

    % Remove cases where patient does not respond
    alphaLSC = [WTalphaLSC resalphaLSC];

    comp_matLSC = [comp_mat(:,1:4) alphaLSC];

    % Filter out nonresponders
    comp_matLSC = comp_matLSC(comp_matLSC(:,5)>pY,:);

    %% Save Results

    csvwrite(string(strcat(drug,'AlphaGeneratorResults121618.csv')),comp_matLSC);

end
    
toc