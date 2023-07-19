

load('psi_fc_data_info','f','info','sub','pla','psi')

% load each of the saved frequency specific corr mats:
freqs = {'DeltaFC' 'ThetaFC' 'AlphaFC' 'BetaFC' 'Gamma1FC' 'Gamma2FC'};
%load DeltaFC.mat

load HCP360.mat

for k = 1:6;
    load(freqs{k});

    % combine placebo and drug for matrix factorisation
    M = [plaRdelta(:,:);psiRdelta(:,:)];
    
    % threshold
    np = 0.1; %-top to keep
    [~,I] = maxpoints(sum(abs(M)),np*(360.^2)); % get indices
    b = find(~ismember(1:(360.^2),I));
    M(:,b) = 0;
    
    % work in reduced (subnetwork) space:
    V  = computereducedoperator(M);
    rM = M*V';    
    K  = 8;
    
    % non-negative matrix factorisation: data-drive identification of
    % "subnetworks"
    [W,H] = nnmf(rM,K);
    
    % mean, group and repeated measures (within sub)
    X = [ones(30,1) [-1*ones(15,1);ones(15,1)] [eye(15);eye(15)] ];

    % weighted least squares estimator - is the minimum-variance unbiased
    % estimator (Gauss-Markov theorem)
    iC = pinv(cov(X));
    iC(1,:) = 1; iC(:,1) = 1;
    b = (pinv(X*iC*X')*X*iC)'*rM;

    % lsq fit
    %b = X\rM;
    G = X*b;
    r = rM - G;
    %R2(k) = 1 - (norm(r)/norm(rM-mean(rM)))^2;
    
    % Produce the subnetworks
    SN = G*H';
    % %[W,H] = nnmf(G,K);

    % find hubs in subnetworks? - i.e. nodes with most connections
    for i = 1:K
        sn{i} = reshape(mean(W(:,i)*H(i,:),1)*V,[360 360]);
    end

    % this routine (definesubnet.m) finds a hub, then maps all components connected to it 
    for i = 1:K
        X0 = sn{i};
        %v = roicentres(reduced.v,reduced.vi);
        [indices{i},subnet{i}] = definesubnet(X0);
    end

    % generate full space subnet masks
    for ik = 1:K
        mask = zeros(360);
        ind = indices{ik};

        for i = 1:size(ind,1)
            x = ind(i,1);
            y = ind(i,2);
            mask(x,y) = 1;
        end
        subnetmask{ik} = mask;

    end

    % reconstrct subnets
    for i = 1:size(W,1)       
        for ik = 1:K

            X = reshape(W(i,ik)*H(ik,:)*V,[360 360]);   
            recon(i,ik,:,:) = X.*subnetmask{ik};

        end
    end

    % assessment of between groups
    R = kRandTest(SN(16:end,:)-SN(1:15,:),[],0.05,5000);
    
    stats{k} = R;

    % generate means for plot
    for i = 1:K
        CompCh{i} = squeeze( mean( recon(16:end,i,:,:) - recon(1:15,i,:,:),1) );
    end

    % plot
    figure('Position',[1225         217        2800        1364]) 
    
    sp = ceil(K/2);
    for i = 1:K
        s = subplot(2,sp,i);
        net = squeeze(CompCh{i});
        atemplate('mesh','def1','sourcemodel',{reduced.v reduced.vi},'network',net,'nocolbar','fighnd',s);
    end

    drawnow;
    export_fig(['MeanChangeDrugMinusPla_' freqs{k} '.png'],'-m3','-transparent')

    
    % % reconstruct topographies of subnetworks
    % for i = 1:K
    % 
    %     C = W(:,i)*H(i,:)*V;
    % 
    %     CP = C(1:15,:);
    %     CD = C(16:end,:);
    % 
    %     for j = 1:15
    % 
    %         % individual network topographies
    %         intopP(j,:,:) = reshape(CP(j,:),[360 360]);
    %         intopD(j,:,:) = reshape(CD(j,:),[360 360]);
    % 
    %     end
    % 
    %     % averages per group for plotting mean change
    %     MeanPlaComp(i,:,:) = squeeze(mean(intopP,1));
    %     MeanPsiComp(i,:,:) = squeeze(mean(intopD,1));
    % 
    % 
    % end

    % try definesubnet
    
    
    % figure('Position',[1225         217        2800        1364])
    % 
    % 
    % sp = ceil(K/2);
    % for i = 1:K
    %     s = subplot(2,sp,i);
    %     net = squeeze(MeanPsiComp(i,:,:)) - squeeze(MeanPlaComp(i,:,:));
    %     net = thresh(net,.1);
    %     atemplate('mesh','def1','sourcemodel',{reduced.v reduced.vi},'network',net,'nocolbar','fighnd',s);
    % end
    % 
    % drawnow;
    % export_fig(['MeanChangeDrugMinusPla_' freqs{k} '.png'],'-m3','-transparent')

end

save('ResultsPSI_BestSoFar')