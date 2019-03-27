%%% Model for transmission bottlenecks in infectious disease


function nres = InfDisease030919(delt)

% clear
% close all

rng(10);
% delt = 3e-5; 

%% Parameters

%%% Rates

% Rates - within host

% Birth rates [/day]
b_S = 0.4;
b_A = b_S;
b_B = b_S;

% Death rates without drug [/day]
d0_S = 0.2;
d0_A = d0_S;
d0_B = d0_S;

% Net growth rates (without drug)
r0_S = b_S-d0_S;
r0_A = b_A-d0_A;
r0_B = b_B-d0_B;

% Drug 1 kill rates [/day]
alpha1_S = 0.3;
alpha1_A = 0.01;
alpha1_B = 0.19;

d1_S = d0_S+alpha1_S;
d1_A = d0_A+alpha1_A;
d1_B = d0_B+alpha1_B;

r1_S = b_S-d1_S;
r1_A = b_A-d1_A;
r1_B = b_B-d1_B;

% Drug 2 kill rates [/day]
alpha2_S = 0.3;
alpha2_A = 0.3;
alpha2_B = 0.3;

d2_S = d0_S+alpha2_S;
d2_A = d0_A+alpha2_A;
d2_B = d0_B+alpha2_B;

r2_S = b_S-d2_S;
r2_A = b_A-d2_A;
r2_B = b_B-d2_B;

% Mutation rates
mut = 1e-8;     % resistance mutation rate [/division]
rho_A = 0.25;    % allele A mutation probability [/resistance mutation]
rho_B = 0.75;    % allele B mutation probability [/resistance mutation]

% Summarize within-host dynamics
pars = [b_S d0_S d1_S d2_S NaN
        b_A d0_A d1_A d2_A rho_A
        b_B d0_B d1_B d2_B rho_B];

% Rates - between hosts
beta = 1e-2;   % contact rate [/day]
% del = 1e-3; % fractional bottleneck size

%%% Populations

N = 500;     % number of hosts

M = 1e7;    % number of pathogens at detection


%% Simulations

nsims = 5e3;

nres = zeros(nsims,3);

Scol = [42 186 252]/255;
Acol = [46 54 143]/255;
Bcol = [236 32 39]/255;

parfor r = 1:nsims
    
    % Describe patient state
    pat_mat = zeros(N,2);
    % col 1 = infection state (0 = S; 1 = I; 2 = R)
    % col 2 = drug treatment state (0 = no drug; 1 = 1st gen; 2 = 2nd gen)

    % Describe pathogenic population state
    pop_mat = zeros(N,3);
    % col 1 = number of pathogens - WT
    % col 2 = number of pathogens - Allele A
    % col 3 = number of pathogens - Allele B

    % Describe time to next within-patient event
    tevt_mat = NaN(N,8);
    % col 1 = time to detection of primary disease
    % col 2 = time to detection of A
    % col 3 = time to detection of B
    % col 4 = time to resistance A
    % col 5 = time to resistance B
    % col 6 = time to elimination of WT (defined as collapse to one cell)
    % col 7 = time to elimination of A
    % col 8 = time to elimination of B
    
    % Patient 1 begins with 100 pathogens
    pat_mat(1,1) = 1;   % update state to infected
    pop_mat(1,1) = 100; % update pathogen population

    % Time of detection for patient 1
    tevt_mat(1,1) = t_det(1,1,0,pars,pop_mat,pat_mat,M);

    % Time of mutation for patient 1 (before detection)
    t_curr = 0;
    idx_pat = 1;
    tevt_mat(1,4) = t_mut(idx_pat,2,[t_curr min(tevt_mat(idx_pat,:))],pars,pop_mat,pat_mat,mut);
    tevt_mat(1,5) = t_mut(idx_pat,3,[t_curr min(tevt_mat(idx_pat,:))],pars,pop_mat,pat_mat,mut);

    % Initialize
    n = 1;
    tn = 0;

    res_det = zeros(N,1);

    trans_order = zeros(1,N);
    trans_order(1) = 1;

    don_node = [];
    acc_node = [];

    idx_G = 1;
%     edgeCol = zeros(0,3);
%     nodeColTot = zeros(0,3);
%     nodeColObs = zeros(0,3);
    
    % Run simulation until no more infected individuals
    while sum(pat_mat(:,1,end)==1)>0

        % Define current state
        pat_curr = pat_mat(:,:,end);
        pop_curr = pop_mat(:,:,end);
        t_curr = sum(tn);    

        % Define net growth rates for next interval
        drug0pats = pat_curr(:,2)==0; % Patients not on treatment
        drug1pats = pat_curr(:,2)==1; % Patients on drug 1
        drug2pats = pat_curr(:,2)==2; % Patients on drug 2

        rate_mat = NaN(N,3);
        rate_mat(drug0pats,:) = repmat([r0_S r0_A r0_B],sum(drug0pats),1); % growth rates for those not on treatment
        rate_mat(drug1pats,:) = repmat([r1_S r1_A r1_B],sum(drug1pats),1); % growth rate for those on drug 1
        rate_mat(drug2pats,:) = repmat([r2_S r2_A r2_B],sum(drug2pats),1); % growth rate for those on drug 2

        %%% Determine next event and inter-event time

        % Transmission events
        Spats = pat_curr(:,1)==0;
        nS = sum(Spats);

        Ipats = pat_curr(:,1)==1;
        nI = sum(Ipats);

        transRate = beta/N*nS*nI;
        timeElapsed_trans = exprnd(1/transRate);

        % Check if transmission event precedes predetermined macroscopic events
        [timeElapsed_det,evt_next] = min(tevt_mat(:)-t_curr);

        if timeElapsed_trans<timeElapsed_det % transmission event precedes predetermined

            % Update time vector
            timeElapsed = timeElapsed_trans;
            t_next = t_curr+timeElapsed;

            % Update population matrix to before event
            pop_next = pop_curr.*exp(rate_mat*timeElapsed);
            pat_next = pat_curr;

            pat_acc = find(Spats,1);

            if nI==1
                pat_don = find(Ipats);
            else
                pat_don = randsample(find(Ipats),1);
            end

            pop_don = pop_next(pat_don,:);

            % Draw founding population according to bottleneck size
            pop_fnd = zeros(1,3);
            for i = 1:3
                pop_i = round(pop_don(i));
                if pop_i*delt < 15
                    pop_fnd(i) = binornd(pop_i,delt);
                else
                    pop_fnd(i) = normrnd(pop_i*delt,sqrt(pop_i*delt*(1-delt)));
                end
            end
            pop_fnd = round(pop_fnd);

            if sum(pop_fnd)>0
                trans_order(idx_G+1)=pat_acc;

                % Update state
                pat_next(pat_acc,1) = 1;        % designate infected
                pop_next(pat_acc,:) = pop_fnd;  % update population

                % Record transmission
                don_node(idx_G) = pat_don;
                acc_node(idx_G) = pat_acc;

                % Sensitive detection of new host
                if pop_fnd(1)>0
                    tevt_mat(pat_acc,1) = t_det(pat_acc,1,t_next,pars,pop_next,pat_next,M);

%                     edgeCol(idx_G,1) = 1;
                end

                % Resistant A dynamics - have to split up if statement to not
                % mess up logic
                if pop_fnd(2)>0 % if A transmitted
                    tevt_mat(pat_acc,2) = t_det(pat_acc,2,t_next,pars,pop_next,pat_next,M); % time of A detection

%                     edgeCol(idx_G,2) = 1;
                end

                % Resistant B dynamics
                if pop_fnd(3)>0 % if B transmitted
                    tevt_mat(pat_acc,3) = t_det(pat_acc,3,t_next,pars,pop_next,pat_next,M); % time of B detection

%                     edgeCol(idx_G,3) = 1;
                else
                    t_nnext = min(tevt_mat(pat_acc,:));
                    tevt_mat(pat_acc,5) = t_mut(pat_acc,3,[t_next t_nnext],pars,pop_next,pat_next,mut); % time of B mutation
                end

                if pop_fnd(2)==0
                    t_nnext = min(tevt_mat(pat_acc,:));
                    tevt_mat(pat_acc,4) = t_mut(pat_acc,2,[t_next t_nnext],pars,pop_next,pat_next,mut); % time of A mutation
                end

%                 if pop_fnd(1)>0&&sum(pop_fnd([2 3]))==0
%                     edgeCol2(idx_G,:) = Scol;
%                 elseif pop_fnd(2)>0&&sum(pop_fnd([1 3]))==0
%                     edgeCol2(idx_G,:) = Acol;
%                 elseif pop_fnd(3)>0&&sum(pop_fnd([1 2]))==0
%                     edgeCol2(idx_G,:) = Bcol;
%                 elseif sum(pop_fnd([1 2]))>0&&pop_fnd(3)==0
%                     edgeCol2(idx_G,:) = mean([Scol;Acol]);
%                 elseif sum(pop_fnd([1 3]))>0&&pop_fnd(2)==0
%                     edgeCol2(idx_G,:) = mean([Scol;Bcol]);
%                 elseif sum(pop_fnd([2 3]))>0&&pop_fnd(1)==0
%                     edgeCol2(idx_G,:) = mean([Acol;Bcol]);
%                 else
%                     edgeCol2(idx_G,:) = [0 0 0];
%                 end

                idx_G = idx_G + 1;

            end

        else

            % Update time vector
            timeElapsed = timeElapsed_det;
            t_next = t_curr+timeElapsed;

            % Update population matrix to before event
            pop_next = pop_curr.*exp(rate_mat*timeElapsed);
            pat_next = pat_curr;

            % Identify event and update state
            [idx_pat,idx_evt] = ind2sub(size(tevt_mat),evt_next);

            % Event: detection of sensitive disease
            if idx_evt==1

                pat_next(idx_pat,2) = 1; % patient begins drug 1            
                tevt_mat(idx_pat,6) = t_ext(idx_pat,1,t_next,pars,pop_next,pat_next); % time of extinction of sensitive population

                if pop_next(idx_pat,2)>0 % if patient has Allele A
                    tevt_mat(idx_pat,2) = t_det(idx_pat,2,t_next,pars,pop_next,pat_next,M); % time of Allele A detection
                else
                    evts = 1:8;
                    t_nnext = min(tevt_mat(idx_pat,evts(evts~=idx_evt)));
                    tevt_mat(idx_pat,4) = t_mut(idx_pat,2,[t_next t_nnext],pars,pop_next,pat_next,mut); % time of Allele A mutation
                end

                if pop_next(idx_pat,3)>0 % if patient has Allele B
                    tevt_mat(idx_pat,3) = t_det(idx_pat,3,t_next,pars,pop_next,pat_next,M); % time of Allele B detection
                else
                    evts = 1:8;
                    t_nnext = min(tevt_mat(idx_pat,evts(evts~=idx_evt)));
                    tevt_mat(idx_pat,5) = t_mut(idx_pat,3,[t_next t_nnext],pars,pop_next,pat_next,mut); % time of Allele B mutation
                end

            % Event: Allele A detected
            elseif idx_evt==2

                pat_next(idx_pat,2) = 2; % patient begins drug 2
                tevt_mat(idx_pat,7) = t_ext(idx_pat,2,t_next,pars,pop_next,pat_next); % time of extinction of A

                if pop_next(idx_pat,1)>0 % if patient has sensitive cells
                    tevt_mat(idx_pat,1) = NaN; % sensitive won't reach detection on drug 2
                    tevt_mat(idx_pat,6) = t_ext(idx_pat,1,t_next,pars,pop_next,pat_next); % update time of extinction of sensitive
                    if pop_next(idx_pat,3)==0 % if patient has no B
                        evts = 1:8;
                        t_nnext = min(tevt_mat(idx_pat,evts(evts~=idx_evt)));
                        tevt_mat(idx_pat,5) = t_mut(idx_pat,3,[t_next t_nnext],pars,pop_next,pat_next,mut); % update time of B mutation
                    end
                end

                if pop_next(idx_pat,3)>0 % if patient has B
                    tevt_mat(idx_pat,3) = NaN; % B won't reach detection on drug 2
                    tevt_mat(idx_pat,8) = t_ext(idx_pat,3,t_next,pars,pop_next,pat_next); % time of B extinction
                end

                res_det(idx_pat) = 2; % note detection of resistance

            % Event: Allele B detected
            elseif idx_evt==3

                pat_next(idx_pat,2) = 2; % patient begins drug 2
                tevt_mat(idx_pat,8) = t_ext(idx_pat,3,t_next,pars,pop_next,pat_next); % time of extinction of B

                if pop_next(idx_pat,1)>0 % if patient has sensitive cells
                    tevt_mat(idx_pat,1) = NaN; % sensitive won't reach detection on drug 2
                    tevt_mat(idx_pat,6) = t_ext(idx_pat,1,t_next,pars,pop_next,pat_next); % update time of extinction of sensitive
                    if pop_next(idx_pat,2)==0 % if patient has no A
                        evts = 1:8;
                        t_nnext = min(tevt_mat(idx_pat,evts(evts~=idx_evt)));
                        tevt_mat(idx_pat,4) = t_mut(idx_pat,2,[t_next t_nnext],pars,pop_next,pat_next,mut); % update time of A mutation
                    end
                end

                if pop_next(idx_pat,2)>0 % if patient has A
                    tevt_mat(idx_pat,2) = NaN; % A won't reach detection on drug 2
                    tevt_mat(idx_pat,7) = t_ext(idx_pat,2,t_next,pars,pop_next,pat_next); % time of A extinction
                end

                res_det(idx_pat) = 3; % note detection of resistance

            % Event: resistance mutation A 
            elseif idx_evt==4

                pop_next(idx_pat,2) = 1; % resistance A spawns
                tevt_mat(idx_pat,2) = t_det(idx_pat,2,t_next,pars,pop_next,pat_next,M); % time of allele A detection

                if pop_next(idx_pat,3)==0 % if patient has no B
                    evts = 1:8;
                    t_nnext = min(tevt_mat(idx_pat,evts(evts~=idx_evt)));
                    tevt_mat(idx_pat,5) = t_mut(idx_pat,3,[t_next t_nnext],pars,pop_next,pat_next,mut); % time of B mutation
                end

            % Event: resistance mutation B
            elseif idx_evt==5

                pop_next(idx_pat,3) = 1; % resistance B spawns
                tevt_mat(idx_pat,3) = t_det(idx_pat,3,t_next,pars,pop_next,pat_next,M); % time of allele B detection

                if pop_next(idx_pat,2)==0 % if patient has no A
                    evts = 1:8;
                    t_nnext = min(tevt_mat(idx_pat,evts(evts~=idx_evt)));
                    tevt_mat(idx_pat,4) = t_mut(idx_pat,2,[t_next t_nnext],pars,pop_next,pat_next,mut); % time of A mutation
                end

            % Event: extinction of any population
            elseif ismember(idx_evt,6:8)

                pop_next(idx_pat,idx_evt-5) = 0; % extinction of population

                if sum(pop_next(idx_pat,:))==0 % if no remaining pathogens
                    pat_next(idx_pat,1) = 2; % patient recovered

                    % Record within-patient population history
                    for i = 1:3
                        if sum(pop_mat(idx_pat,i,:))>0 % if patient had variant
%                             nodeColTot(idx_pat,i) = 1;
                        end
                    end
                end

%                 if idx_evt==7
%                     nodeColObs(idx_pat,:) = Acol;
%                 elseif idx_evt==8
%                     nodeColObs(idx_pat,:) = Bcol;
%                 end

            end

            tevt_mat(idx_pat,idx_evt) = NaN;    % remove evt from event-time mat

        end

        % Update patient matrix

        % Store results
        tn(n+1) = timeElapsed;
        pat_mat(:,:,n+1) = pat_next;
        pop_mat(:,:,n+1) = pop_next;

        n = n+1;

    end

    nR = sum(pat_mat(:,1,end)==2);

    nA = sum(res_det==2);
    nB = sum(res_det==3);
    nS = nR - nA - nB;

    nres(r,:) = [nS nA nB];
    
    disp(r);
    
end

pres = mean(nres);
stdres = std(nres);

%% Plot results

% % for i = 1:length(find(trans_order))
% for i = 1:4
%     figure(i)
%     hold on
%     for j = 1:3
%         plot(cumsum(tn),reshape(pop_mat(trans_order(i),j,:),1,[]))
%     end
%     hold off
% end

% G = digraph(don_node,acc_node);
% 
% figure
% p = plot(G,'Layout','force');

% % replace white nodes with black
% whiteRows = find(sum(nodeColTot,2)==3);
% if whiteRows>0
%     nodeColTot(whiteRows,:) = repmat([0 0 0],length(whiteRows),1);
% end
% p.NodeColor = nodeColTot;
% 
% % reorder edgeCol
% p.EdgeColor = edgeCol(G.Edges.EndNodes(:,2)-1,:);

% figure
% p2 = plot(G,'Layout','force');
% 
% % Node Colors
% noResRows = find(sum(nodeColObs,2)==0);
% nodeColObs(noResRows,:) = repmat(Scol,length(noResRows),1);
% nodeColObs = [nodeColObs;
%               repmat(Scol,size(G.Nodes,1)-size(nodeColObs,1),1)];
% p2.NodeColor = nodeColObs;
% 
% % Edge Colors
% p2.EdgeColor = edgeCol2(G.Edges.EndNodes(:,2)-1,:);
% 
% p2.NodeLabel={};