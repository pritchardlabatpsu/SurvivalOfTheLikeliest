%%% Spatial model of tumor growth

clear
close all

rng(15);

%%% TO DO:
% spatial parameters

%% Parameters

%%% Time Parameters

% Each ball = 1e5 cells --> scale rates appropriately
sf = 1e6;

% Rates
b = 0.69;       % [/day] birth rate (Waclaw 2015)
d = 0.35;       % [/day] death rate (Waclaw 2015)
c = 0.5;        % [/day] clearance rate (~24 hour removal - Elliott Dev Cell 2017)
mut = 1e-8*sf;  % [/division] resistance mutation rate

% Drug kill rates
netgrowth_post = d-b; % (Waclaw 2015)
alphatreatS = b-d-netgrowth_post;
alphatreatA = 0.15*(b-d);
alphatreatB = 0.85*(b-d);

% Mutation probabilities
rhoA = 0.2;    % mutation bias A [/resistance mutation]
rhoB = 0.8;    % mutation bias B [/resistance mutation]

%%% Spatial parameters

% Spatial gradient constants
kb = -10;       % Diffusion constant (birth)
kd = 10;        % Diffusion constant (death)
kalpha = -10;   % Diffusion constant (drug killing)

%%% Simulation Parameters
nsims = 1e3;
res_profile = NaN(nsims,4);

tic

parfor z = 1:nsims

    %%% Populations
    
    M = 5e8/sf;
    cellmax = 2*M;
    
    pos_curr = NaN(cellmax,3);
    stt_curr = NaN(cellmax,2);

    % Position matrix
    pos_curr(1,:,1) = [0 0 0]; % [x,y,z]

    % State matrix
    % col 1: alive = 1, dead = 0
    % col 2: genotype - S = 1, A = 2, B = 3
    stt_curr(1,:) = [1 1];

    liv = [];
    ded = [];

    %% Simulation

    nghbrs1 = combvec(-1:1,-1:1,-1:1)';
    [~,row0s] = ismember([0 0 0],nghbrs1,'rows');
    nghbrs1(row0s,:) = [];

    % Initialize
    n = 1;
    t = 0;

    ncells_next = 1;
    nviab_next = 1;

    treat = false;
    response = false;
    alphaS = 0;
    alphaA = 0;
    alphaB = 0;

    tic
    % Run simulation until:
    % no remaining viable cells or
    % relapse after response or
    % no response (treatment failure)
    while nviab_next>0 && (~response || ncells_next<M) && ncells_next<cellmax

        % Update current state
%         pos_curr = pos_mat(:,:,end);
%         stt_curr = stt_mat(:,:,end);

        cells = pos_curr(~isnan(pos_curr(:,1)),:);
        state = stt_curr(~isnan(pos_curr(:,1)),:);
        ncells = size(cells,1);

        viab = cells(state(:,1)==1,:);
        nviab = size(viab,1);

        dead = cells(state(:,1)==0,:);
        ndead = size(dead,1);

        viabS = cells(ismember(state,[1 1],'rows'),:);
        nviabS = size(viabS,1);

        viabA = cells(ismember(state,[1 2],'rows'),:);
        nviabA = size(viabA,1);

        viabB = cells(ismember(state,[1 3],'rows'),:);
        nviabB = size(viabB,1);

        pos_next = pos_curr;
        stt_next = stt_curr;

        % Treatment state
        if ncells>=M
            treat=true;
            nM = n;
        end

        if treat
            alphaS = alphatreatS;
            alphaA = alphatreatA;
            alphaB = alphatreatB;
        end

        % Determine event
        evts = [b*nviab         % division
                d*nviab         % natural death event
                alphaS*nviabS   % drug kill event (S)
                alphaA*nviabA   % drug kill event (A)
                alphaB*nviabB   % drug kill event (B)
                c*ndead];       % drug clearance event
        theta = sum(evts);
        t(n) = -log(rand)/theta;                            % time to next event
        idx_evt = randsample(length(evts), 1, true, evts);  % event

        if idx_evt==1 % division event

            %%% Determine cell i to divide
            
            % Identify distance of cells from periphery
            if ncells>30
                bndry = boundary(cells,1);
                [vrtx_idx,~,rnk] = unique(bndry(:));
                vrtcs = cells(vrtx_idx,:);
                bndry_rnk = reshape(rnk,size(bndry));

                % Find distance of each live cell from tumor edge
                prphdst = point2trimesh('Faces',bndry_rnk,'Vertices',vrtcs, 'QueryPoints', viab); 

            else
                prphdst = zeros(nviab,1);
            end

            % Find distance
            idxi = randsample(nviab,1,true,exp(kb*(abs(prphdst))));
            posi = viab(idxi,:);
            stti = state(ismember(cells,posi,'rows'),:);

            %%% Determine relative position of new cell j

            % Define "cubic sphere" with radius k around cell i and 
            % scan for free spaces
            posf = NaN;
            k = 1;
            prevscan = [0 0 0];

            while isnan(posf)

                % Begin with cube of length 2k
                nghbrhdk = combvec(-k:k,-k:k,-k:k)';

                % Remove previously scanned positions and those more than k
                % units from position i
                nghbrhdk(ismember(nghbrhdk,prevscan,'rows'),:) = [];

                normsk = sqrt(nghbrhdk(:,1).^2+nghbrhdk(:,2).^2+nghbrhdk(:,3).^2);
                nghbrhdk(normsk>k,:) = [];

                % Scan for positions k units from position i
                nghbrsik = posi+nghbrhdk;
                freek = nghbrsik(~ismember(nghbrsik,cells,'rows'),:);

                freeik = freek-posi;

                % Find nearest free space (position f) to cell i
                normsik = sqrt(freeik(:,1).^2+freeik(:,2).^2+freeik(:,3).^2);
                minfreek = find(normsik==min(normsik));

                if isempty(minfreek)
                    prevscan = [prevscan; nghbrhdk];
                    k = k+1;
                elseif length(minfreek)==1
                    posf = freek(minfreek,:);
                elseif length(minfreek)>1
                    posf = freek(randsample(minfreek,1),:);
                end

            end

            % Push cells between i and f until free space adjacent to i
            % Definitions:
            % position i = site of dividing cell
            % position f = site of initially free space
            % position g = site of cell that moves to position f
            % position j = site of new cell (adjacent to position i)

            nghbrsi1 = posi + nghbrs1;

            while ~ismember(posf,nghbrsi1,'rows')

                % Find position g adjacent to f closest to line through i
                % and f

                vec_if = posf - posi;
                dstr = NaN(26,1);
                for r = 1:26
                    % Candidate for migrant cell g
                    posgr = posf+nghbrs1(r,:);
                    vec_ig = posgr-posi;

                    % Consider only positions g that are closer to i than f
                    % is to i
                    if norm(vec_ig) < norm(vec_if)
                        % Use parallelogram definition of cross product to
                        % find distance between candidate position g and
                        % line passing through i and f
                        dstr(r) = norm(cross(vec_ig,vec_if))/norm(vec_if);
                    end
                end

                % Choose position g
                idx_mindst = find(dstr==min(dstr));
                if length(idx_mindst)>1
                    posg = posf+nghbrs1(randsample(idx_mindst,1),:);
                else
                    posg = posf+nghbrs1(idx_mindst,:);
                end

                % Move cell in position g to position f
                idx_posg = ismember(pos_next,posg,'rows');
                pos_next(idx_posg,:) = posf;

                % Position g is now open to receive next pushed cell
                posf = posg;

            end

            % Once a free space is adjacent to cell i, it divides ==> cell j
            posj = posf;

            % Determine genotype of cell j
            if rand<mut && stti(2)==1
                if rand<rhoA
                    sttj = [1 2];
                else
                    sttj = [1 3];
                end
            else
                sttj = stti;
            end

            % Update by filling in NaN rows
            idxj = find(all(isnan(pos_next),2),1,'first');
            
            % May be necessary to expand matrices
            if isempty(idxj)
                idxj = size(pos_next,1)+1;
            end
            
            pos_next(idxj,:) = posj;
            stt_next(idxj,:) = sttj;

        elseif idx_evt==2 % natural death event

            %%% Determine cell i to die
            
            % Identify distance of cells from periphery
            if ncells>30
                bndry = boundary(cells,1);
                [vrtx_idx,~,rnk] = unique(bndry(:));
                vrtcs = cells(vrtx_idx,:);
                bndry_rnk = reshape(rnk,size(bndry));

                % Find distance of each live cell from tumor edge
                prphdst = point2trimesh('Faces',bndry_rnk,'Vertices',vrtcs, 'QueryPoints', viab); 

            else
                prphdst = zeros(nviab,1);
            end
            
            % Choose cell i
            idxi = randsample(nviab,1,true,exp(kd*abs(prphdst)));
            posi = viab(idxi,:);

            % Update
            stt_next(ismember(pos_next,posi,'rows'),1) = 0;

        elseif idx_evt==3 % drug death of S

            %%% Determine cell i to die
            
            % Identify distance of cells from periphery
            if ncells>30
                bndry = boundary(cells,1);
                [vrtx_idx,~,rnk] = unique(bndry(:));
                vrtcs = cells(vrtx_idx,:);
                bndry_rnk = reshape(rnk,size(bndry));

                % Find distance of each live cell from tumor edge
                prphdstS = point2trimesh('Faces',bndry_rnk,'Vertices',vrtcs, 'QueryPoints', viabS); 

            else
                prphdstS = zeros(nviabS,1);
            end
            
            % Choose cell i
            idxi = randsample(nviabS,1,true,exp(kalpha*(abs(prphdstS))));
            posi = viabS(idxi,:);

            % Update
            stt_next(ismember(pos_next,posi,'rows'),1) = 0;

        elseif idx_evt==4 % drug death of A
            
            % Identify distance of cells from periphery
            if ncells>30
                bndry = boundary(cells,1);
                [vrtx_idx,~,rnk] = unique(bndry(:));
                vrtcs = cells(vrtx_idx,:);
                bndry_rnk = reshape(rnk,size(bndry));

                % Find distance of each live cell from tumor edge
                prphdstA = point2trimesh('Faces',bndry_rnk,'Vertices',vrtcs, 'QueryPoints', viabA); 

            else
                prphdstA = zeros(nviabA,1);
            end
            
            % Determine cell i to die
            idxi = randsample(nviabA,1,true,exp(kalpha*(abs(prphdstA))));
            posi = viabA(idxi,:);

            % Update
            stt_next(ismember(pos_next,posi,'rows'),1) = 0;

        elseif idx_evt==5 % drug death of B

            %%% Determine cell i to die
            
            % Identify distance of cells from periphery
            if ncells>30
                bndry = boundary(cells,1);
                [vrtx_idx,~,rnk] = unique(bndry(:));
                vrtcs = cells(vrtx_idx,:);
                bndry_rnk = reshape(rnk,size(bndry));

                % Find distance of each live cell from tumor edge
                prphdstB = point2trimesh('Faces',bndry_rnk,'Vertices',vrtcs, 'QueryPoints', viabB); 

            else
                prphdstB = zeros(nviabB,1);
            end
            
            % Choose cell i
            idxi = randsample(nviabB,1,true,exp(kalpha*(abs(prphdstB))));
            posi = viabB(idxi,:);

            % Update
            stt_next(ismember(pos_next,posi,'rows'),1) = 0;

        elseif idx_evt==6 % clearance of dead cell

            % Determine cell i to be removed
            idxi = randsample(ndead,1);
            posi = dead(idxi,:);

            % Update
            stt_next(ismember(pos_next,posi,'rows'),:) = NaN(1,2);
            pos_next(ismember(pos_next,posi,'rows'),:) = NaN(1,3);

            %%% Push cells from periphery to fill in empty space

            % Define "cubic sphere" with radius k around position i and 
            % scan for free spaces
            posf = NaN;
            k = 1;
            prevscan = [0 0 0];

            while isnan(posf)

                % Begin with cube of length 2k
                nghbrhdk = combvec(-k:k,-k:k,-k:k)';

                % Remove previously scanned positions and those more than k
                % units from position i
                nghbrhdk(ismember(nghbrhdk,prevscan,'rows'),:) = [];

                normsk = sqrt(nghbrhdk(:,1).^2+nghbrhdk(:,2).^2+nghbrhdk(:,3).^2);
                nghbrhdk(normsk>k,:) = [];

                % Scan for positions k units from position i
                nghbrsik = posi+nghbrhdk;
                freek = nghbrsik(~ismember(nghbrsik,cells,'rows'),:);

                freeik = freek-posi;

                % Find nearest free space (position f) to cell i
                normsik = sqrt(freeik(:,1).^2+freeik(:,2).^2+freeik(:,3).^2);
                minfreek = find(normsik==min(normsik));

                if isempty(minfreek)
                    prevscan = [prevscan; nghbrhdk];
                    k = k+1;
                elseif length(minfreek)==1
                    posf = freek(minfreek,:);
                elseif length(minfreek)>1
                    posf = freek(randsample(minfreek,1),:);
                end

            end

    %         %%% ALT
    %         
    %         % Find vertices of tumor
    %         bndry = boundary(cells,1);
    %         vrtx_idx = unique(bndry(:));
    %         vrtcs = cells(vrtx_idx,:);
    %         
    %         % Identify vertex closest to removed cell
    %         posiv = vrtcs-posi;
    %         normsiv = sqrt(posiv(:,1).^2+posiv(:,2).^2+posiv(:,3).^2);
    %         
    %         mindstv = find(normsiv==min(normsiv));
    %         if length(mindstv)==1
    %             posf = posi+posiv(mindstv);
    %         else
    %             posf = posi+posiv(randsample(mindstv,1),:);
    %         end
    %         
            %%%

            % Push cells between f and i
            % Definitions:
            % position i = site of internal free space (removed cell)
            % position f = site of peripheral free space
            % position g = site of cell that moves to position i

            nghbrsf1 = posf + nghbrs1;
            nghbrsf0 = [nghbrsf1; posf];

            while ~ismember(posi,nghbrsf0,'rows')

                % Find position g adjacent to i closest to line through i
                % and f

                vec_fi = posi - posf;
                dstr = NaN(26,1);
                for r = 1:26
                    % Candidate for migrant cell g
                    posgr = posi+nghbrs1(r,:);
                    vec_fg = posgr-posf;

                    % Consider only positions g that are closer to f than i
                    % is to f
                    if norm(vec_fg) < norm(vec_fi)
                        % Use parallelogram definition of cross product to
                        % find distance between candidate position g and
                        % line passing through i and f
                        dstr(r) = norm(cross(vec_fg,vec_fi))/norm(vec_fi);
                    end
                end

                % Choose position g
                idx_mindst = find(dstr==min(dstr));
                if length(idx_mindst)>1
                    posg = posi+nghbrs1(randsample(idx_mindst,1),:);
                else
                    posg = posi+nghbrs1(idx_mindst,:);
                end

                % Move cell in position g to position i
                idx_posg = ismember(pos_next,posg,'rows');
                pos_next(idx_posg,:) = posi;

                % Position g is now open to receive next pushed cell
                posi = posg;

            end

        end

        % Update matrices
%         pos_mat(:,:,n+1) = pos_next;
%         stt_mat(:,:,n+1) = stt_next;
        
        pos_curr = pos_next;
        stt_curr = stt_next;

%         S(n) = nviabS;
%         A(n) = nviabA;
%         B(n) = nviabB;
% 
%         liv(n) = nviab;
%         ded(n) = ndead;

        ncells_next = sum(~isnan(pos_next(:,1)));
        nviab_next = sum(stt_next(:,1)==1);
%         disp(ncells_next)
        
        % Note when patient responds
        if treat && ncells_next<=0.9*M
            response = true;
        end

        n = n+1;

    end
        
    % Record results if tumor reached detection
    if treat
        res_profile(z,:) = [nviabS nviabA nviabB ndead];
    end
    disp(z)
    
end

toc

%% Save Results

filename = 'SpatialModelLoopResults032419.csv';
csvwrite(filename,res_profile);

%% Plotting

% % Plot limits
% minx = min(min(pos_mat(:,1,:)));
% maxx = max(max(pos_mat(:,1,:)));
% miny = min(min(pos_mat(:,2,:)));
% maxy = max(max(pos_mat(:,2,:)));
% minz = min(min(pos_mat(:,3,:)));
% maxz = max(max(pos_mat(:,3,:)));
% 
% % Colors
% % Alive
% Sacol = [42 186 252]/255;
% Aacol = [46 54 143]/255;
% Bacol = [236 32 39]/255;
% 
% % Dead
% Sdcol = [156 216 244]/255;
% Adcol = [93 97 137]/255;
% Bdcol = [229 144 147]/255;
% 
% clrs = NaN(cellmax,3);
% 
% figure
% for i = 1:n
%     
%     % Colors
%     Sarows = ismember(stt_mat(:,:,i),[1 1],'rows');
%     Aarows = ismember(stt_mat(:,:,i),[1 2],'rows');
%     Barows = ismember(stt_mat(:,:,i),[1 3],'rows');
%     
%     Sdrows = ismember(stt_mat(:,:,i),[0 1],'rows');
%     Adrows = ismember(stt_mat(:,:,i),[0 2],'rows');
%     Bdrows = ismember(stt_mat(:,:,i),[0 3],'rows');
%     
%     clrs(Sarows,:) = repmat(Sacol,sum(Sarows),1);
%     clrs(Aarows,:) = repmat(Aacol,sum(Aarows),1);
%     clrs(Barows,:) = repmat(Bacol,sum(Barows),1);
%     
%     clrs(Sdrows,:) = repmat(Sdcol,sum(Sdrows),1);
%     clrs(Adrows,:) = repmat(Adcol,sum(Adrows),1);
%     clrs(Bdrows,:) = repmat(Bdcol,sum(Bdrows),1);
%         
%     if i==nM
%         clrsM = clrs;
%     end
%     
%     scatter3(pos_mat(:,1,i),pos_mat(:,2,i),pos_mat(:,3,i),1e3,clrs,'.');
%     
%     xlim([minx-1 maxx+1]);
%     ylim([miny-1 maxy+1]);
%     zlim([minz-1 maxz+1]);
%     view(-60,25)
% %     set(gca,'Color','k')
%     pause(1e-8)
% end
% 
% time = cumsum(t);
% 
% figure
% plot(time/365,liv*sf)
% hold on
% plot(time/365,ded*sf)

% figure
% P = pos_mat(:,:,nM);
% xc = P(P(:,3)==0,:);
% clrs2d = clrsM(P(:,3)==0,:);
% scatter3(xc(:,1),xc(:,2),xc(:,3),1e3,clrs2d,'.')
% view(2)
% 
% figure
% scatter3(P(:,1),P(:,2),P(:,3),1e3,clrsM,'.');
% view(-60,25)
% 
% k = boundary(P,1);
% hold on
% trisurf(k,P(:,1),P(:,2),P(:,3),'Facecolor','blue','FaceAlpha',0.1);
% 
% faces = boundary(P,1);
% [vert_idx,~,rnk] = unique(faces(:));
% vertices = P(vert_idx,:);
% face_rnk = reshape(rnk,size(faces));
% points = viab;
% 
% figure
% hold on
% [distances,surface_points] = point2trimesh('Faces',face_rnk,'Vertices',vertices, 'QueryPoints', points); 
% %  patch('Faces',faces','Vertices',vertices,'FaceAlpha',.5); xlabel('x'); ylabel('y'); zlabel('z'); axis equal; hold on
% trisurf(face_rnk,vertices(:,1),vertices(:,2),vertices(:,3),'Facecolor','red','FaceAlpha',0.1);
% plot3M = @(XYZ,varargin) plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),varargin{:});
% plot3M(points,'*r')
% plot3M(surface_points,'*k')
% plot3M(reshape([shiftdim(points,-1);shiftdim(surface_points,-1);shiftdim(points,-1)*NaN],[],3),'k')
% view(-60,25)
% 
