
function output = ATSims022719(pop_treat,pop_remain,pars)

% clear
% close all

rng(25);

%%% Inputs

% pop_treat = 1e9;
% pop_remain = 1e6;
% 
% pars = csvread('ATParameters022419.csv');
biasR = pars(:,1);
rel_alpha = pars(:,2);

%%% Constant Parameters %%%

% Populations
WT_0 = 100;    % initial sensitive population [cells]

Rsites = length(biasR);
N = Rsites + 1; % number of variant populations to be modeled (resistant variants + WT)
pop_init = zeros(1,N);
pop_init(1) = WT_0;

% Rates [/day] (Bozic 2013)
b = 0.14;
d = 0.13;
alpha_WT = 0.04; % may decrease to a more conservative estimate

alpha_res = (b-d)*rel_alpha;

mut = 1e-8;

%%% Model Parameters %%%

tau_init = NaN;     % length of model time intervals [days] (if NaN model determines tau adaptively)
epsilon = 0.03;     % error tolerance []

npats = 5e4;

%% Simulations

% Initialize storage variables

% timeSet = cell(npats,1);
% cellSet = cell(npats,1);

output = zeros(npats,1);

parfor r = 1:npats

    % Initialize
    n = 1;
    t = 0;
    tau = tau_init;
    singleSSAs = NaN;
    treat = false;
    cells = pop_init;
    alpha = zeros(1,N);
    alpha_treat = [alpha_WT alpha_res'];

    % Run simulation until WBC population reaches detection (relapse)
    while sum(cells(end,:))>0 && sum(cells(end,2:end))<pop_treat

        % Define current
        P_curr = cells(end,:);
        P_WT = P_curr(1);
        P_res = P_curr(2:end);
        t_curr = sum(t);

        % Treatment conditions
        if sum(P_curr) >= pop_treat && ~treat
            treat = true;
            
            % Debulk tumor
            P_next = zeros(1,N);
            P_next(1) = pop_remain;
            P_next(2:end) = binornd(P_curr(2:end),pop_remain/pop_treat);
            
            cells(n+1,:) = P_next;
            
            timeElapsed = 0;
            t(n+1) = timeElapsed;
            
            n = n+1;
            continue         
        end

        if treat
            alpha = alpha_treat;
        end      

        % Define propensity vector and change matrix  
        evts = [b*(1-mut)*P_WT           % Faithful WT division
                b*P_res'                % Resistant division
                ((d+alpha).*P_curr)'    % Cell death
                b*mut*P_WT*biasR];       % Mutation
        
        numChange = [1 zeros(1,Rsites)              % Faithful WT division
                     zeros(Rsites,1) eye(Rsites)    % Resistant division
                     -eye(N)                        % Cell death
                     zeros(Rsites,1) eye(Rsites)];  % Mutation
        
        % Calculate total rate
        theta = sum(evts);

        % Tau-leaping model
        if isnan(singleSSAs)

            % Determine tau if not defined
            if isnan(tau)

                % Evaluate Jacobian matrix at P_curr
                Jac = [b*(1-mut) zeros(1,Rsites)         % Faithful WT division
                       zeros(Rsites,1) b*eye(Rsites)    % Resistant division
                       diag(d+alpha)                    % Cell death
                       b*mut*biasR zeros(Rsites)];       % Mutation

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
    
    if sum(cells(end,:))>0
        [~,domAllele] = max(cells(end,2:end));
        output(r) = domAllele;
    end
    
%     timeSet{r} = cumsum(t);
%     cellSet{r} = cells;


    disp(r);
end

% % Save results
% csvwrite('ATSims022419.csv',output);

%% Plot Results

% colors = [0 0 1; rand(Rsites,3)];
% 
% figure
% for i = 1:npats
%     for j = 1:N
%         semilogy(timeSet{i}/365,cellSet{i}(:,j),'Color',colors(j,:))
%         hold on
%     end
% end

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