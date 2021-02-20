function [Wexcerpt , regIDX , sCorMat , meanFR] = dynamic_clustering_morpho(WM , muSTD , somConst , outString)
    T = 25920;%000/3; % full morpology only for 5 days, otherwise too computationally expensive
    % import tree morphology with TREES toolbox
    resConst = 10; %resampling constant (micron)
    t = load_tree('../data/L23.swc');
    % resample tree into equally spaced segments
    tr = resample_tree(t , resConst);
    % get matrix of all daugther segments
    dauMat = getDaughterMat(tr);
    % only basal and apical dendrite
    subTree = find(tr.R == 4 | tr.R == 3); 
    den = 0.2;
    N = floor(den*length(subTree)*resConst);
    % get adjacency matrix and compute shortest path distances between all
    % segments
    sAdj = tr.dA + tr.dA';
    denDisMat = distances(graph(sAdj))*resConst;
    % model parameters
    beta = 0.5; taum = 600;  eta = 0.45;  tauy = 300;
    Winit = 0.5; tauW = 3000/50 * (1/(2*beta*(1 - eta))); 
    tauY = tauy/50; tauR = taum/50; 
    rho = (2*eta - 1)/(2 - 2*eta); 
    phi = 3;eps = 0.02; normSTD = 6; simID = now;
    % bAP parameters
    bAP_Thresh = 25; bAP_AMP = 5;
    % LN model parameters
    a = 0.01; b = 9.4;
    % distribute synapses by selecting random segments and random positions
    % within segments
    pos = subTree(randsample(1:length(subTree)  ,  N  , 1));
    subpos = -rand(N , 1)*resConst;
    % the actual distance is then computed as the distance between segments
    % plus or minus the distances within the segments
    dMat = denDisMat(pos , pos);
    dMat = augDMat(dMat , dauMat(pos , pos) , subpos);
    % transform into proximity matrix S
    SMat = normpdf(dMat , 0 , normSTD)*(sqrt(2*pi)*normSTD);
    SMat = SMat - SMat.*eye(N) + eye(N);
    % compute attenuation factor for each segment
    compSomDist = PL_tree(tr)*resConst;
    somDVec = normpdf(compSomDist(pos) , 0 , somConst)*(sqrt(2*pi)*somConst);
    % randomly sample receptive field centers. Factor 0.8 to ensure that
    % the entire RF is within visual space
    RFDIST = trandn(-0.8*(pi/muSTD)*ones(N,1),0.8*(pi/muSTD)*ones(N,1))*muSTD; % first sample the RF offset
    RFANGLE = rand(N,1)*2*pi - pi; % then the angle in visual space
    [MUx , MUy] = pol2cart(RFANGLE , RFDIST);
    MUs = [MUx , MUy];
    % sample orientation preferences
    thetas = rand(N,1)*2*pi;
    % generate 2D Gabors
    fullVisualFieldDegree = 100;
    S1 = 1*(2*pi)/fullVisualFieldDegree; 
    S2 = 2*(2*pi)/fullVisualFieldDegree; % receptive field sizes in custom coordinate system
    M = size(WM , 1);
    [XSPACE , YSPACE] = meshgrid(linspace(-pi , pi , M) , linspace(-pi , pi , M));
    fxs = zeros( M , M , N);
    for ii = 1:N
        [fx] = mixGauss([XSPACE(:) , YSPACE(:)] , MUs(ii , :) , thetas(ii) , S1 , S2);
        fxs(: , : , ii) = reshape(fx , [M , M]);
    end
    % set up accumulators
    R = zeros(N,T); Y = zeros(N,T); YnoBAP = zeros(N,T);
    W = zeros(N,T); S = zeros(N,T); W(:,1) = Winit;
    A = zeros(1 , T);  Ain = 0; regIDX = [];
    storeConfig = struct(); storeConfig.thetas = [thetas]; storeConfig.pos = [pos]; storeConfig.W = [W(:,1)]; storeConfig.MUs = MUs;
    % iterate through time
    for tt = 2:T
        if mod(tt , 10000) == 0; fprintf('%d,',tt); end
        if mod(tt , 100000) == 0; fprintf('\n'); end
        % get current retinal wave frame
        cWM = WM(: , : , mod(ceil((tt-1)/2) , size(WM , 3) - 1) + 1);
        % generate input events
        [Sin , ~] = getSin(cWM , fxs , a , b);  S(: , tt) = Sin;

        R(: , tt) = R(: , tt-1)*exp(-1./tauR) + phi*Sin*(1 - exp(-1./tauR));
        Waug = repmat(W(: , tt-1)' , [N , 1]).*SMat; 
        Y(: , tt) = Y(: , tt-1)*exp(-1./tauY) + (Waug*Sin + bAP_AMP*(Ain).*somDVec)*(1 - exp(-1./tauY)); 
        % also compute postsynaptic accumulator without bAP
        YnoBAP(: , tt) = YnoBAP(: , tt-1) +  -(1./tauY)*(YnoBAP(: , tt-1) -  Waug*Sin);
        W(: , tt) =  min(max(W(: , tt-1) + (1./tauW)*( Y(: , tt-1).*  (R(: , tt-1) + rho) ) , 0),1) ; % compute weight dynamics
        % compute the somatic activation from the synaptic accumulators
        % without bAP
        A( tt) = sum(W(: , tt).*YnoBAP(: , tt)); Ain = (A(tt-1) > bAP_Thresh)*(rand < 1/4);
        % if synaptic efficacy below cutoff
        if sum(W(: , tt) < eps) > 0
            % resample
            pos(W(:,tt) < eps) = subTree(randsample(1:length(subTree)  ,  1  , 1));
            subpos(W(:,tt) < eps) = -rand(1 , 1)*resConst;
            dMat = denDisMat(pos , pos);
            dMat = augDMat(dMat , dauMat(pos , pos) , subpos);
            SMat = normpdf(dMat , 0 , normSTD)*(sqrt(2*pi)*normSTD);
            SMat = SMat - SMat.*eye(N) + eye(N);
            somDVec = normpdf(compSomDist(pos) , 0 , somConst)*(sqrt(2*pi)*somConst);
            replID = find(W(: , tt) < eps);
            for kk = 1:length(replID)
                RFDIST = trandn(-0.8*(pi/muSTD),0.8*(pi/muSTD))*muSTD;
                RFANGLE = rand(1,1)*2*pi - pi;
                [MUx , MUy] = pol2cart(RFANGLE , RFDIST);
                MUs(replID(kk) , 1) = MUx; MUs(replID(kk) , 2) = MUy;
                thetas(replID(kk)) = rand(1,1)*2*pi;
                [fx] = mixGauss([XSPACE(:) , YSPACE(:)] , MUs(replID(kk) , :) , thetas(replID(kk)) , S1 , S2);
                fxs(: , : , replID(kk)) = reshape(fx , [M , M]);   
                regIDX = [regIDX ; tt];
            end
            storeConfig.thetas = [ storeConfig.thetas , thetas]; storeConfig.pos = [ storeConfig.pos , pos];
            storeConfig.W = [storeConfig.W , W(:,tt)]; storeConfig.MUs = cat(3 , storeConfig.MUs , MUs);
            W(W(:,tt) < eps , tt-1) = nan; W(W(:,tt) < eps , tt) = 0.5;
        end
    end
    fprintf('\n')
    Wexcerpt = W(: , end-10000:end); Rexcerpt = R(: , end-10000:end); 
    Yexcerpt = Y(: , end-10000:end); Sexcerpt = S(: , end-100000:end);
    Aexcerpt = A(end-10000:end);
    [~ , sID] = sort(pos); sCorMat = corrcoef(smoothdata(S(sID , tt-10000:tt)' , 'movmean' ,  20)); meanFR = mean(S(sID , tt-10000:tt) , 2);
    clear S R W Y WM A YnoBAP
    save(sprintf('../sims/%s_%0.20f.mat' , outString , simID) )
end
