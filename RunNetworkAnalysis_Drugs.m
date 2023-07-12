% This script will load up the precomputed functional connectivity matrices
% for each dataset and frequency band and run the main analyses.
%
% The aim is to look at functional connectivity networks that are common
% across the ket/pla and psi/pla datasets, and compare them statistically.
% Sounds simple enough...
%
% There are 4 main steps:
% 1  - initial thresholding top 20% of connections across all datasets
% 2a - use NNMF to derive subnetworks (needs choice of N, I'm using 8)
% 2b - reduce each of these K subnetworks to a set of only interconnected
%      network nodes (i.e. a network consisting a hub and all connected
%      nodes but all other removed)
% 3  - generate design matrix with mean, drug/no drug, ket only, psi only
%      and within-subjects columns.
% 4  - use W-LSQ to fit model to subnetworks (not raw data) and assess 
%      statistics (fit, predictors). Note this is ONE GLM incorporating all
%      subnetworks.
% 5  - project individual data through subnetworks and run randomisation 
%      based t-tests on each subnetwork separately.
% 6  - VISUALISE mean effects
%
% AS2023

% note this script requires 'aconnectivity' toolbox for analysis functions
% and atemplate (SourceMesh toolbox) for plotting.

method = 'nnmf'; % can be nnmf, svd, eig or pca, qr


cd('/Users/alexandershaw/Library/CloudStorage/Dropbox/PSI_KET_2023/');

% load each of the saved frequency specific corr mats:
freqs = {'DeltaFC' 'ThetaFC' 'AlphaFC' 'BetaFC' 'Gamma1FC' 'Gamma2FC'};

load HCP360.mat
v = aconnectivity.roicentres(reduced.v,reduced.vi);

for k = 1:6;
    PSI = load(['PSILO/' freqs{k}]);
    KET = load(['KET/' freqs{k}]);

    % combine placebo and drug for matrix factorisation
    %--------------------------------------------------------------------
    f0 = fieldnames(PSI);
    f1 = fieldnames(KET);
    M = [PSI.(f0{contains(f0,'pla')})(:,:); PSI.(f0{contains(f0,'psi')})(:,:);
         KET.(f1{contains(f1,'pla')})(:,:); KET.(f1{contains(f1,'ket')})(:,:) ];
    
    % indices for groups / studies:
    psipla = 1:15;
    psipsi = 16:30;
    ketpla = 31:47;
    ketket = 48:64;
    
    % STEP 1: initial threshold
    %--------------------------------------------------------------------
    np = 0.1; %-top to keep/2
    [~,I] = aconnectivity.maxpoints(sum(abs(M)),np*(360.^2)); % get indices
    b = find(~ismember(1:(360.^2),I));
    M(:,b) = 0;
    
    % work in reduced (thresh-network) space: rM = M*V
    V  = aconnectivity.computereducedoperator(M,@var);
    rM = M*V';   

    K  = 6;
    
    % STEP 2: non-negative matrix factorisation: 
    % data-driven identification of "subnetworks"
    %--------------------------------------------------------------------
    if K > 1
        switch method
            case 'nnmf';
                [W,H] = nnmf(rM,K);%,'replicates',20);
            case 'svd';
                [W,S,H] = svd(rM);
    
                W = W(:,1:K);
                S = S(1:K,1:K);
                H = H(:,1:K)';
                W = W*S;
    
            case 'eig'
                [H,W] = eigs(cov(rM),K);
                H = H';
    
            case 'pca';
                [W,S,H] = svd(cov(rM),'econ');
    
                W = W(:,1:K);
                S = S(1:K,1:K);
                H = H(:,1:K)';
                W = W*S;
            case 'qr'
                [q,r] = qr(rM);

                W = q(:,1:K);
                H = r(1:K,:);

            case 'orth'
                H = orth(rM')';
                H = H(1:K,:);

            case 'pod';
                [W,S,H] = pod(rM');
        
                S = diag(S);
                W = W(:,1:K);
                S = S(1:K,1:K);
                H = H(:,1:K)';
                W = W*S;
            case 'lu'
                [L,U] = lu(rM);
                W = L(:,1:K);
                H = U(1:K,:);
                
        end
    else
        H = ones(size(rM,2),1)';
    end

    % reconstruct subnetworks in full space
    for i = 1:K
        sn{i} = reshape(H(i,:)*V,[360 360]);
        sn{i} = (sn{i} + sn{i}')./2;
    end

    % map all components connected to a hub
    %--------------------------------------------------------------------
    for i = 1:K
        X0 = sn{i};   

        [indices{i},subnet{i},clustersize(i)] = aconnectivity.definesubnet(X0,[],[],.98);
        subnet{i} = (subnet{i}+subnet{i}')./2;

        % if missed, re-try with no theresholding
        if clustersize(i) == 0
            [indices{i},subnet{i},clustersize(i)] = aconnectivity.definesubnet(X0,[],[],1);
            subnet{i} = (subnet{i}+subnet{i}')./2;
        end
    end

    % reverse create right side of nnmf from subnetworks
    for i = 1:K
        subnet{i}  = (subnet{i}+subnet{i}')/2;
        HH(i,:)    = V*subnet{i}(:);
        REGLABS{i} = labels(find(sum(reshape(HH(i,:)*V,[360 360]))));
    end
    
    % normalise subnetwork weights
    %HH = HH./norm(HH);

    % STEP 3: set up design matrix for regression/GLM/ANOVA model
    %--------------------------------------------------------------------
    X = [ones(64,1) [0*ones(15,1);ones(15,1);0*ones(17,1);ones(17,1)] ...
        [zeros(15,1);ones(15,1);zeros(17,1);zeros(17,1)]...
        [zeros(15,1);zeros(15,1);zeros(17,1);ones(17,1)]...
        [[eye(15);eye(15)];zeros(17*2,15)] [zeros(15*2,17);eye(17);eye(17)] ];
    

    % STEP 4: Fit the model using weighted least squares estimator
    % the minimum-variance unbiased estimator (Gauss-Markov theorem)
    %--------------------------------------------------------------------
    [b,G,stats] = aconnectivity.aglm(X,rM*HH');

    % NOTE: this "full model" is in subnetwork space, so we are testing for
    % predictor effects not on all connections but K-different networks

    % GLM ASSESSMENT:
    %-------------------------------------------------------------------
    % was the GLM fit significant?
    fprintf('Full GLM Model Fit: F = %d, p = %d\n',stats.F,stats.p);

    % is there an effect of drug?
    fprintf('Was there an effect of drug? F = %d, p = %d\n',stats.Xtval(2),stats.Xpval(2));

    % is there an effect of psilo?
    fprintf('Was there an effect of psilo? F = %d, p = %d\n',stats.Xtval(3),stats.Xpval(3));
    
    % is there an effect of ketamine?
    fprintf('Was there an effect of ketamine? F = %d, p = %d\n',stats.Xtval(4),stats.Xpval(4));


    % assessment of drugs on different subnetworks - i.e. post-hoc tests
    %--------------------------------------------------------------------
    % PSI - PLA
    Rpsi = aconnectivity.kRandTest(G(16:30,:)-G(1:15,:),[],0.05,5000);

    % KET - PLA
    Rket = aconnectivity.kRandTest(G(48:64,:)-G(31:47,:),[],0.05,5000);

    % RMAN: (PSI-PLA) vs (KET-PLA)
    Rdrug = aconnectivity.kRandTest(G(16:30,:)-G(1:15,:),G(48:64,:)-G(31:47,:),0.05,5000);

    % save statistics
    astats(k).Rpsi  = Rpsi;
    astats(k).Rket  = Rket;
    astats(k).Rdrug = Rdrug;

    % save also matrix factorisation and subnetwork
    astats(k).W     = W;
    astats(k).H     = H;
    astats(k).HH    = HH;
    
    % save also GLM stat data
    astats(k).stats = stats;
    astats(k).subjcomps = rM*HH';
    astats(k).X = X;
    astats(k).full = X*b*HH;


    % compute relative subnet contribution for each subject
    %projection = X*b*HH;
    %for is = 1:size(rM,1)
    %    bs(:,is) = HH'\projection(is,:)';
    %end
    %[hh,pp,~,sta] = ttest(bs(:,drug)',bs(:,pla)');

    % display all the subnets by colour on a brain
    everynet = zeros(360);
    w = [-3 -2 -1 1 2 3];
    for i = 1:K
        Z = w(i) * ~~subnet{i};
        everynet = everynet + Z;
    end
    afigure,atemplate('mesh','def1','sourcemodel',{reduced.v reduced.vi},...
        'network',everynet,'nodes',sum(everynet));

    export_fig(['AllSubnets_' freqs{k} '.png'],'-m3','-transparent');
    drawnow;
    close;


    % second analysis: how does drug change subnetwork weighting? in terms
    % of probability: i.e. what is the probability of being in each
    % subnetwork for each dataset?
    %----------------------------------------------------------------------
    P = (W' ./ sum(W'))';

    [bP,PP,Pstats] = aconnectivity.aglm(X,P);

    % assessment of drugs on different subnetworks - i.e. post-hoc tests
    %--------------------------------------------------------------------
    % PSI - PLA
    Ppsi = aconnectivity.kRandTest(PP(16:30,:)-PP(1:15,:),[],0.05,5000);

    % KET - PLA
    Pket = aconnectivity.kRandTest(PP(48:64,:)-PP(31:47,:),[],0.05,5000);

    % RMAN: (PSI-PLA) vs (KET-PLA)
    Pdrug = aconnectivity.kRandTest(PP(16:30,:)-PP(1:15,:),PP(48:64,:)-PP(31:47,:),0.05,5000);

    % save statistics
    astats(k).ProbabSubnet_si  = Ppsi;
    astats(k).ProbabSubnet_ket  = Pket;
    astats(k).ProbabSubnet_drug = Pdrug;

    % create plots
    figure,hold on;

    plot(1:6,mean(PP(psipla,:),1),'linewidth',3);
    plot(1:6,mean(PP(psipsi,:),1),'linewidth',3);
    plot(1:6,mean(PP(ketpla,:),1),'linewidth',3);
    plot(1:6,mean(PP(ketket,:),1),'linewidth',3);

    legend({'pPLA' 'PSI' 'kPLA' 'KET'});

    export_fig(['ProbabilityOfSubnetByGroup_' freqs{k} '.png'],'-m3','-transparent');


    % generate individual networks and drug means for plot
    %--------------------------------------------------------------------
    for i = 1:K

        % the subject-by-region-by-region map
        netmap{i} = reshape( W(:,i)*HH(i,:)*V, [size(G,1) 360 360]);
        
        % generate mean change values for plots
        CompChKpsi{i} = squeeze( mean( netmap{i}(16:30,:,:) - netmap{i}(1:15,:,:),1) );
        CompChKket{i} = squeeze( mean( netmap{i}(48:64,:,:) - netmap{i}(31:47,:,:),1) );

    end

    % plot
    %--------------------------------------------------------------------
    figure('Position',[1225         217        2800        1364]) 
    fname = {'Delta' 'Theta' 'Alpha' 'Beta' 'Gamma1' 'Gamma2'};

    sp = K;%ceil(K/2);
    for i = 1:K
        s = subplot(3,sp,i);
        
        % PSI plot
        net0 = squeeze(CompChKpsi{i});
        atemplate('mesh','def1','sourcemodel',{reduced.v reduced.vi},'network',net0,'fighnd',s,'nodes',sum(net0),'netcmap',alexcmap);
        
        t = astats(k).Rpsi.tseries_corr(i);
        p = astats(k).Rpsi.pseries_corr(i);

        title(sprintf([fname{k} ' C%d PSI v PLA\nt=%d, p=%d'],i,t,p));
        set(findall(gca, 'type', 'text'), 'visible', 'on');

        % KET plot
        s = subplot(3,sp,i+K);
        net1 = squeeze(CompChKket{i});
        atemplate('mesh','def1','sourcemodel',{reduced.v reduced.vi},'network',net1,'fighnd',s,'nodes',sum(net0),'netcmap',alexcmap);

        t = astats(k).Rket.tseries_corr(i);
        p = astats(k).Rket.pseries_corr(i);

        title(sprintf([fname{k} ' C%d KET v PLA\nt=%d, p=%d'],i,t,p));
        set(findall(gca, 'type', 'text'), 'visible', 'on');

        % DIFF PLOT
        s = subplot(3,sp,i+(2*K));
        netx = net0 - net1;
        atemplate('mesh','def1','sourcemodel',{reduced.v reduced.vi},'network',netx,'fighnd',s,'nodes',sum(net0),'netcmap',alexcmap);

        t = astats(k).Rdrug.tseries_corr(i);
        p = astats(k).Rdrug.pseries_corr(i);

        title(sprintf([fname{k} ' C%d PSI v KET\nt=%d, p=%d'],i,t,p));
        set(findall(gca, 'type', 'text'), 'visible', 'on');

    end

    drawnow;
    export_fig(['MeanChangeDrugeffects_' freqs{k} '.png'],'-m3','-transparent')


end

save('ResultsPSI_BestSoFar')