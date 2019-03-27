%%%% Mutational Bias in Predicting Resistance Outgrowth
%%%% Scott Leighow 1/26/18

function output_res = LSCMutBiasGeneral_SL011619(alpha_mat,biasR)

rng(25);

%%% Inputs

pars = alpha_mat(:,1:3);
alphas = alpha_mat(:,4:end);

%%% Constant Parameters %%%

% Populations
Y0_WT = 50;    % initial sensitive LSC population [cells]

Rsites = length(biasR);

M = Rsites + 1; % number of variant populations to be modeled (resistant variants + WT)
N = 3*M;        % total number of subpopulations to be modeled

%           [ QLSCs  ] [      PLSCs     ] [  WBCs  ]
pop_init = [zeros(1,M) Y0_WT zeros(1,M-1) zeros(1,M)];  % initial population structure [cells]

% Rates [/day] (Fassoni 2018)
pXY = 0.05380/30.4; % QLSC activation rate [/day]
pYX = 0.00095/30.4; % PLSC deactivation rate [/day]
pY = 0.18701/30.4;  % PLSC net growth (without drug) [/day]
q = 1.00394/30.4;   % PLSC net drug kill term [/day]
pW = 4e6/30.4;      % PLSC differentiation [/day]
rW = 4/30.4;        % WBC death [/day]

% Break down PLSC net growth (pY) into birth and death rates
turnover = 0.3;     % estimate from Komarova and Wodarz 2005
bY = pY/(1-turnover);
dY = bY-pY;

mut = 4e-8; % Michor 2005

%%% Model Parameters %%%

tau_init = NaN;     % length of model time intervals [days] (if NaN model determines tau adaptively)
epsilon = 0.03;     % error tolerance []
niter = 1;

%% Simulations

% Initialize storage variables
% outMat = zeros(npats,5+4*Rsites);

% timeSet = cell(npats,1);
% cellSet = cell(npats,1);

nrow = size(alpha_mat,1);

output = NaN(nrow,niter,5);

tic

parfor r = 1:nrow
    
    output_r = [repmat(pars(r,:),niter,1) (1:niter)' NaN(niter,1)];
    
    for iter = 1:niter

        % Initialize
        n = 1;
        t = 0;
        tau = tau_init;
        singleSSAs = NaN;
        treat = false;
        cells = pop_init;
        alpha = zeros(1,M);
        alpha_treat = alphas(r,:);
        
        pop_treat = 10^pars(r,3);

        % Run simulation until WBC population reaches detection (relapse)
        while sum(cells(end,[M+1 N-Rsites+1:N]))>0 && sum(cells(end,N-Rsites+1:N))<pop_treat


            % Define current
            P_curr = cells(end,:);        
            t_curr = sum(t);

            X = P_curr(1:M);
            Y = P_curr(M+1:2*M);
            W = P_curr(2*M+1:3*M);

            % Treatment conditions
            if sum(P_curr) >= pop_treat
                treat = true;
            end

            if treat
                alpha = alpha_treat;
            end      

            % Define propensity vector and change matrix        
            evts = [pXY*X'               % QLSC activation
                    pYX*Y'               % PLSC deactivation
                    bY*(1-mut)*Y(1)      % faithful WT PLSC division
                    bY*Y(2:end)'         % resistant PLSC division
                    ((dY+alpha).*Y)'     % PLSC death
                    pW*Y'                % PLSC differentiation
                    rW*W'                % WBC death
                    bY*mut*Y(1)*biasR']; % resistance mutation

            numChange = [-eye(M)      eye(M)                  zeros(M)      % QLSC activation
                         eye(M)       -eye(M)                 zeros(M)      % PLSC deactivation
                         zeros(M)     eye(M)                  zeros(M)      % PLSC division
                         zeros(M)     -eye(M)                 zeros(M)      % PLSC death
                         zeros(M)     zeros(M)                eye(M)        % PLSC differentiation
                         zeros(M)     zeros(M)                -eye(M)       % WBC death
                         zeros(M-1,M) zeros(M-1,1) eye(M-1)   zeros(M-1,M)];% resistance mutation

            % Calculate total rate
            theta = sum(evts);

            % Tau-leaping model
            if isnan(singleSSAs)

                % Determine tau if not defined
                if isnan(tau)

                    % Evaluate Jacobian matrix at P_curr
                    Jac = [pXY*eye(M)   zeros(M)                 zeros(M)        % QLSC activation
                           zeros(M)     pYX*eye(M)               zeros(M)        % PLSC deactivation
                           zeros(1,M)   bY*(1-mut) zeros(1,M-1)  zeros(1,M)      % faithful WT PLSC division
                           zeros(M-1,M) zeros(M-1,1) bY*eye(M-1) zeros(M-1,M)    % resistant PLSC division
                           zeros(M)     diag(dY+alpha)           zeros(M)        % PLSC death
                           zeros(M)     pW*eye(M)                zeros(M)        % PLSC differentiation
                           zeros(M)     zeros(M)                 rW*eye(M)       % WBC death
                           zeros(M-1,M) bY*mut*biasR' zeros(M-1) zeros(M-1,M)];  % resistance mutation


                    % Calculate tau
                    fjj = Jac*numChange';
                    mu = fjj*evts;
                    std = (fjj.^2)*evts;

                    tau = min(min(epsilon*theta./(abs(mu)), epsilon^2*theta^2./(std)));

                end

                % Rerun with singleSSAs if tau is too small
                if tau < 1/theta
                    singleSSAs = 10;
                    tau = tau_init;
                    continue

                % Otherwise calculate new state
                else
                    % For each event, find number of times it's occured
                    % in time interval tau
                    P_next = P_curr;
                    for i = 1:length(evts)
                        numEvt = poissrnd(evts(i)*tau);
                        P_next = P_next + numChange(i,:).*numEvt;
                    end
                    P_next = round(P_next);

                    % Rerun with smaller tau if any populations are
                    % negative
                    if sum(P_next < 0) > 0
                        tau = tau/2;
                        continue

                    % Otherwise save results in array
                    else                      
                        timeElapsed = tau;
                        t(n+1) = timeElapsed;
                        cells(n+1,:) = P_next;

                        % Reset tau
                        tau = tau_init;
                    end
                end

            % Gillespie without tau leaping
            else
                r1 = rand;
                timeElapsed = -log(r1)/theta;
                t_curr = sum(t(1:n));

                t(n+1) = timeElapsed;

                idx = randsample(length(evts), 1, true, evts);
                P_next = round(P_curr + numChange(idx,:));
                cells(n+1,:) = P_next;

                % Reduce singleSSAs
                singleSSAs = singleSSAs - 1;
                if singleSSAs < 1
                    singleSSAs = NaN;
                end
            end

            n = n + 1;
        end

        % In the case of eradication, repeat if occured pretreatment
        if sum(cells(end,:)) == 0 && ~treat
            continue
        end

        % Save time vector as cumulative sum of time intervals
%         time = cumsum(t)';

    %     figure
    %     semilogy(time,cells)
    %     disp(cells(end,:))

    %     % Sort alleles
    %     relng_wodrug = (bR-d)/(bS-d);
    %     relng_wdrug = (bR-d-alpha1R)/(bS-d);
    %     prop = cells(end,:)/sum(cells(end,:))*100;
    %     
    %     BiasFitProp = [biasR; relng_wodrug; relng_wdrug; prop(2:end)]';
    %     
    %     sortedBiasFitProp = sortrows(BiasFitProp,[-1,-2,-3,-4]);
    %     
    %     % Save results
    %     outvec = reshape(sortedBiasFitProp',1,[]);
    %     
    %     outMat(r,:) = [log10(pop_treat) log10(mut) gene r prop(1) outvec];

    %%%%
    
%         % Save results
%         time = cumsum(t);
%         timeSet{r} = time;
%         cellSet{r} = cells;
%         
%         %%%%

        if sum(cells(end,N-M+2:N))>0
            [~,maxidx] = max(cells(end,N-M+2:N));
            output_r(iter,5) = maxidx;
        else
            output_r(iter,5) = 0;
        end
    
    end
    
    output(r,:,:) = output_r;
       
    disp(r);
end

% pres_obs=sum(domAllele~=0)/npats;

output_res = reshape(output,[],5);

toc

%% Plot Results

% colors = [0 0 1; rand(Rsites,3)];
% 
% figure
% % subplot(1,2,1);
% % for i = 1:npats
% %     for j = 1:M
% %         LSCs = cellSet{i}(:,1:M);
% %         WBCs = cellSet{i}(:,3*M+1:4*M);
% %         semilogy(timeSet{i}/365,LSCs(:,j),'o','Color',colors(j,:))
% %         hold on
% %     end
% % end
% % 
% % title('LSCs')
% % xlabel('Time [years]')
% % ylabel('Population Size [cells]')
% % hold off
% % 
% % subplot(1,2,2);
% for i = 1:npats
%     for j = 1:M
%         LSCs = cellSet{i}(:,1:M);
%         WBCs = cellSet{i}(:,3*M+1:4*M);
%         semilogy(timeSet{i}/365,LSCs(:,j),'o','Color',colors(j,:))
%         hold on
%     end
% end
% title('White Blood Cell Dynamics')
% xlabel('Time [years]')
% ylabel('Population Size [cells]')
% hold off
% 
% figure
% for j = 1:Rsites
%     plot(biasR(j),IC50R(j),'.','Color',colors(j+1,:),'markers',25)
%     hold on
% end
% 
% title('Parameter Space')
% xlabel('Probability')
% ylabel('IC_{50} [nM]')
% hold off
% 
% max(cellSet{i}(:,1))
