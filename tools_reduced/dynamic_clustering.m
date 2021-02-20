function [Wexcerpt , regIDX , sCorMat , meanFR] = dynamic_clustering(WM , muSTD , outString)
    %Model constants
    T = 25920000; 
    beta = 0.5; taum = 600;  eta = 0.45;  tauy = 300;
    Winit = 0.5; tauW = 3000/50 * (1/(2*beta*(1 - eta))); 
    tauY = tauy/50; tauR = taum/50; 
    rho = (2*eta - 1)/(2 - 2*eta);  phi = 3;eps = 0.02; L = 150; 
    simID = now;
    
    % density of synapses, length of branch -> number of synapses
    synDensity = 0.2; L = 150; N = floor(L*synDensity); 
    normSTD = 6; % calcium spread constant

    % LN model parameters
    a = 0.01; b = 9.4;
    % initial position of synapses
    pos = sort(rand(N,1))*L; sPos = mod((pos + L/2),L);
    dMat = min(pdist2(pos,pos) , pdist2(sPos,sPos));
    dMat = normpdf(dMat , 0 , normSTD)*(sqrt(2*pi)*normSTD);
    % initial receptive field centers
    RFDIST = trandn(-0.8*(pi/muSTD)*ones(N,1),0.8*(pi/muSTD)*ones(N,1))*muSTD;
    RFANGLE = rand(N,1)*2*pi - pi;
    [MUx , MUy] = pol2cart(RFANGLE , RFDIST);
    MUs = [MUx , MUy];
    % intial receptive field orientations
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
    % initalize accumulators
    R = zeros(N,T); Y = zeros(N,T); YnoBAP = zeros(N,T);  W = zeros(N,T); S = zeros(N,T);

    W(:,1) = Winit;
    
    % set up storage structure
    storeConfig = struct();
    storeConfig.thetas = [thetas];
    storeConfig.pos = [pos];
    storeConfig.W = [W(:,1)];

    for tt = 2:T
        % logging
        if mod(tt , 10000) == 0; fprintf('%d,',tt); end
        if mod(tt , 100000) == 0; fprintf('\n'); end
        % get current RW input
        cWM = WM(: , : , mod(ceil((tt-1)/2) , size(WM , 3) - 1) + 1);
        % LN model
        [Sin , ~] = getSin(cWM , fxs , a , b); S(: , tt) = Sin;
        % neighborhood weights
        Waug = repmat(W(: , tt-1)' , [N , 1]).*dMat;
        % update variables
        R(: , tt) = R(: , tt-1)*exp(-1./tauR) +...
                    phi*Sin*(1 - exp(-1./tauR)); 
        Y(: , tt) = Y(: , tt-1)*exp(-1./tauY) + (Waug*Sin)*(1 - exp(-1./tauY));
        
        W(: , tt) =  min(max(W(: , tt-1) + (1./tauW)*( Y(: , tt-1).*  (R(: , tt-1) + rho) ) , 0),1);
        
        numTO = sum(W(: , tt) < eps);
        if  numTO > 0
            % generate new position
            pos(W(:,tt) < eps) = rand(numTO,1)*L; sPos = mod((pos + L/2),L);
            dMat = min(pdist2(pos,pos) , pdist2(sPos,sPos));
            dMat = normpdf(dMat , 0 , normSTD)*(sqrt(2*pi)*normSTD);
            % new receptive fields
            replID = find(W(: , tt) < eps);
            for kk = 1:length(replID)
                RFDIST = trandn(-0.8*(pi/muSTD),0.8*(pi/muSTD))*muSTD;
                RFANGLE = rand(1,1)*2*pi - pi;
                [MUx , MUy] = pol2cart(RFANGLE , RFDIST);
                MUs(replID(kk) , 1) = MUx;
                MUs(replID(kk) , 2) = MUy;
                
                thetas(replID(kk)) = rand(1,1)*2*pi;
                [fx] = mixGauss([XSPACE(:) , YSPACE(:)] , MUs(replID(kk) , :) , thetas(replID(kk)) , S1 , S2);
                fxs(: , : , replID(kk)) = reshape(fx , [M , M]);        
            end
            W(W(:,tt) < eps , tt-1) = nan; W(W(:,tt) < eps , tt) = 0.5;
            % store configuration
            storeConfig.thetas = [ storeConfig.thetas , thetas]; storeConfig.pos = [ storeConfig.pos , pos]; storeConfig.W = [storeConfig.W , W(:,tt)];            
        end
    end
    fprintf('\n')
    % get excerpts for storage
    Wexcerpt = W(: , end-10000:end); Rexcerpt = R(: , end-10000:end); Yexcerpt = Y(: , end-10000:end); Sexcerpt = S(: , end-100000:end);
    regIDX = find(isnan(sum(W,1)))*50/(24*60*60*1000); regIDXraw = find(isnan(W));
    [~ , sID] = sort(pos); sCorMat = corrcoef(smoothdata(S(sID , tt-10000:tt)' , 'movmean' ,  20)); meanFR = mean(S(sID , tt-10000:tt) , 2);
    % clear large variables
    clear S R W Y WM 
    save(sprintf('../sims/%s_%0.20f.mat' , outString , simID) )
end

