%% FAS_titrateNumPCs_ampDecode.m
%
%
%

%% perform population decoding

nIter = 200; % Number of model iterations to run
targetFreq = 2; 
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
usePCA = true; % use PCA step prior to model training\
useFreqEffNeu = true; % use neurons with significant anova freq effects
useAllNeu = false; % use all available neurons for modeling (no random cell draws)
cellType = ["PV","SOM","EXC"]; % order for plotting
nCellDraw = 150; % 150 (AmpDecode) and 160 (freqDecode) datset2
nCellType = length(cellType);

% Create output structure for storing PCA information
PC = struct();
for nc = 1 : length(cellType)
    PC.(cellType(nc)) = struct();
end

nPCused = 1:1:20;

% CREATE OUTPUT ARRAYS
% Performance of Cell Types (all)
% PCTall = []; %for populating with accuracy measurements
PCTall = zeros(nIter,nCellType,length(nPCused)); % (nj,ni,nd)
% Models (all)
MDLall = cell(nIter,nCellType,length(nPCused)); % (nj,ni,nd)

for nc = 1:length(cellType) % 
    disp(strcat(" ",cellType(nc)))
    
    trainTrial = round(minTrialNum*trainPortion);
    testTrial = round(minTrialNum*(1-trainPortion));
    trainAmp = repmat(ampValAll,trainTrial,1);
    testAmp = repmat(ampValAll,testTrial,1);
    predAmpAll = zeros(nIter,nCellType,nAmp*testTrial);
    testAmpAll = zeros(nIter,nCellType,nAmp*testTrial);

    if useFreqEffNeu == true
        if multiAmp == true
            idx = find(Tc.identity == cellType(nc) & (is.AmpMod | is.DualMod | is.InterMod));
        end
    elseif useFreqEffNeu == false
        idx = find(Tc.identity == cellType(nc));
    end
    
    % for np = 1 : length(nPCused)
    %     fprintf(strcat(" ",num2str(nPCused(np))));
        
    % reiterate random draws
    for ni = 1:nIter

        % progress tracker
        if rem(ni,nIter/10) == 0
            fprintf(" . ");
        end
        
        % randomly draw cells without replacement
        if useAllNeu == true
            selectedCells = idx;
            nCellDraw = length(idx);
        elseif useAllNeu == false
            selectedCells = randsample(idx, nCellDraw, false);
        end
       
        % pre-allocate training and testing data array
        trainResp = zeros(nFreq,trainTrial,nCellDraw);
        testResp = zeros(nFreq,testTrial,nCellDraw);

        % loop over selected cells
        for nk = 1:length(selectedCells)
            cellResp = Tc.responses{selectedCells(nk)};
            cellFreq = Tc.frequencies{selectedCells(nk)};
            cellAmp = Tc.amplitudes{selectedCells(nk)};
                      
            for na = 1:nAmp
                % randomly draw trials at each amplitude, no replacement
                trialNum = find(cellFreq == freqVal(targetFreq) & cellAmp == ampValAll(na));
                selectedTrials = randsample(trialNum, trainTrial+testTrial, false);

                % populate training and testing arrays
                trainResp(na,:,nk) = cellResp(selectedTrials(1:trainTrial));
                testResp(na,:,nk) = cellResp(selectedTrials(trainTrial+1:end));
            end
        end

        % Reshape matricies
        trainResp = reshape(trainResp, nAmp*trainTrial,nCellDraw);
        testResp = reshape(testResp, nAmp*testTrial,nCellDraw);

        % Z-Scoring 
        for nz = 1 : nCellDraw
            mu = mean(trainResp(:,nz));
            sigma = std(trainResp(:,nz));
            trainResp(:,nz) = (trainResp(:,nz) - mu) ./ sigma;
            testResp(:,nz) = (testResp(:,nz) - mu) ./ sigma;
        end

        % Training WITH PCA 
        [coeff,score,latent] = pca(trainResp);
        PC.(cellType(nc)).coeff{ni} = coeff;
        PC.(cellType(nc)).score{ni} = score;
        PC.(cellType(nc)).latent{ni} = latent;
        % totalVariance = sum(latent);
        % CVA = cumsum((latent./totalVariance).*100); %cumulative variance explained

        for np = 1 : length(nPCused)
            % kComp = floor(trainTrial*nFreq*0.5);
            trainResp_red = score(:,nPCused(1):nPCused(np));
            mdl = fitcdiscr(trainResp_red,trainAmp,'DiscrimType','linear');
            predAmp = predict(mdl,testResp*coeff(:,nPCused(1):nPCused(np)));
    
            % collect output
            PCTall(ni,nc,np) = mean(predAmp == testAmp); %accuracy
            predAmpAll(ni,nc,:) = predAmp;
            testAmpAll(ni,nc,:) = testAmp;
            MDLall{ni,nc,np} = mdl;
        end
    end 

    % end %N PCS USED
    fprintf(newline)
end





%%  PLOT ACCURACY TOGETHER (NEW WAY?)
figure('Color',[1 1 1],'Theme','Light'); hold on;
% plot fake data for correct legend order
plot([-1 -1],[-1 -1],'k-'); 
% plot([-1 -1],[-1 -1],'k--');
fill([-1 -1 -1], [-1 -1 -1],colors.PV,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.SOM,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.EXC,'FaceAlpha',0.5,'EdgeColor','none');
for nt = 1 : size(PCTall,4)
    for nc = 1 : length(cellType)
        yMeans = mean(squeeze(PCTall(:,nc,:,nt)),1).*100;
        % ySems = sem(squeeze(PCTall(:,nc,:,nt)),1).*100; % standard error
        yCIs = zeros(length(nPCused),2); % confidence intervals (95%)
        for np = 1 : length(nPCused)
            yCIs(np,:) = mkCI(squeeze(PCTall(:,nc,np,nt))).*100;
        end
        xLooped = [nPCused flip(nPCused)]; 
        % yLooped = [yMeans+ySems flip(yMeans-ySems)]; when using SEM
        yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
        fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.5,'EdgeColor','none');
        switch nt 
            case 1
                plot(nPCused,yMeans,'-','LineWidth',1,'Color',colors.(cellType(nc)));
            case 2
                plot(nPCused,yMeans,'--','LineWidth',1,'Color',colors.(cellType(nc)));
        end
    end
end

% plot([15 500],[20 20],'r--'); % chance performance level 
set(gca,'XLim',[0.5 40.5],'YLim',[50 75],'YTick',10:10:100);
ylabel('model accuracy (%)','FontName','Arial');
xlabel('PC number','FontName','Arial');
title('Frequency decoding performance by cell type','FontName','Arial');
subtitle('From 280 neurons with frequency selectivity','FontName','Arial');
% legend({"Freq. Mod.","All neurons","PV","SOM","EXC"},'Location','NorthWest','box','off');
legend({"Freq. Mod.","PV","SOM","EXC"},'Location','NorthWest','box','off');




%%  PLOT ACCURACY TOGETHER AS DELTA VALUES

figure('Color',[1 1 1],'Theme','Light'); hold on;
% plot fake data for correct legend order
plot([-1 -1],[-1 -1],'k-'); 
% plot([-1 -1],[-1 -1],'k--');
fill([-1 -1 -1], [-1 -1 -1],colors.PV,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.SOM,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.EXC,'FaceAlpha',0.5,'EdgeColor','none');
for nt = 1 : 1 %size(PCTall,4)
    for nc = 1 : length(cellType)
        yMeans = mean(squeeze(PCTall(:,nc,:,nt)),1).*100;
        yDelta = [yMeans(1) diff(yMeans)];
        % ySems = sem(squeeze(PCTall(:,nc,:,nt)),1).*100; % standard error
        yCIs = zeros(length(nPCused),2); % confidence intervals (95%)
        for np = 1 : length(nPCused)
            yCIs(np,:) = mkCI(squeeze(PCTall(:,nc,np,nt))).*100;
        end
        yCIsFix = yCIs - repmat(yMeans',1,2);
        yCIsDelta = repmat(yDelta',1,2) + yCIsFix;
        xLooped = [nPCused flip(nPCused)]; 
        % yLooped = [yMeans+ySems flip(yMeans-ySems)]; when using SEM
        % yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
        yLooped = [yCIsDelta(:,1)' flip(yCIsDelta(:,2))'];
        fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.5,'EdgeColor','none');
        switch nt 
            case 1
                % plot(nPCused,yMeans,'-','LineWidth',1,'Color',colors.(cellType(nc)));
                plot(nPCused,yDelta,'-','LineWidth',1,'Color',colors.(cellType(nc)));
            case 2
                % plot(nPCused,yMeans,'--','LineWidth',1,'Color',colors.(cellType(nc)));
                plot(nPCused,yDelta,'--','LineWidth',1,'Color',colors.(cellType(nc)));
        end
    end
end

% plot([15 500],[20 20],'r--'); % chance performance level 
% set(gca,'XLim',[1 40],'YLim',[60 100],'YTick',60:10:100);
% set(gca,'XLim',[1 6],'XTick',1:6,'YLim',[-5 80],'YTick',0:10:80);
% ylabel('model accuracy (%)','FontName','Arial');
ylabel('\Delta model accuracy (%)','FontName','Arial');
xlabel('PC number','FontName','Arial');
title('Amplitude decoding performance gain per PC','FontName','Arial');
subtitle('From 150 neurons with amplitude selectivity','FontName','Arial');
% legend({"Freq. Mod.","All neurons","PV","SOM","EXC"},'Location','NorthWest');
% legend({"Freq. Mod.","PV","SOM","EXC"},'Location','NorthWest');
set(gca,'XLim',[1 10],'XTick',1:10,'YLim',[-5 60],'YTick',0:10:80,'XScale','log','Box','off');
% set(gca,'XLim',[2 20],'XTick',[2:2:20],'YLim',[-5 40],'YTick',0:10:80,'XScale','log','Box','off');


%% PLOT INSET FIGURE SHOWING TOP 3 PCs OF EACH CELL TYPE (DATASET 1)

targetRep = 10;  % target model iteration to visualize
pcVisNum = 3; % number of PCs to visualize
figure('theme','light','color',[1 1 1]);
tiledlayout(2,2);
nexttile(1)
spx = [-10 0 10];
ntx = [1 2 4];

for nc = 1 : length(cellType)
    nexttile(ntx(nc)); hold on

    mdl = MDLall{targetRep,nc,pcVisNum}; 
    PCresps = zeros(trainTrial,5,pcVisNum);
    for np = 1 : numPCs
        for na = 1 : nFreq
            PCresps(:,na,np) = mdl.X(mdl.Y == freqVal(na),np);
            yMean = mean(PCresps(:,na,np));
        end
    end

    for np = 1 : pcVisNum
        xVals = freqVal + spx(np); % frequency + xspacing  
        locResps = squeeze(PCresps(:,:,np));
        for na = 1 : nFreq
            plot([xVals(na) xVals(na)],mkCI(locResps(:,na)),'-','Color',colors.(cellType(nc))); % CI bar
            plot(xVals(na),mean(locResps(:,na)),'.','Color',colors.(cellType(nc))) % mean point
        end 
        switch np
            case 1
                plot(xVals,mean(locResps,1),'-','Color',colors.(cellType(nc)));
            case 2
                plot(xVals,mean(locResps,1),'--','Color',colors.(cellType(nc)));
            case 3
                plot(xVals,mean(locResps,1),':','Color',colors.(cellType(nc)));
        end
    end
    % set(gca,'YLim',[-4 5],'YTick',-4:2:4,'XTick',freqVal,'XLim',[280 1120]);
    set(gca,'YLim',[-10 8],'YTick',-8:4:8,'XTick',freqVal,'XLim',[280 1120]);
    ylabel("\DeltaF/F",'FontName','Arial');
    xlabel("Frequency (Hz)",'FontName','Arial');
    hold off
end

% LEGEND IN SPOT 3
nexttile(3); hold on;
for nc = 1 : 3
plot(0,0,'.','Color',colors.(cellType(nc)),'MarkerSize',12);
end
plot([0 0],[0 0],'k-'); 
plot([0 0],[0 0],'k--');
plot([0 0],[0 0],'k:');
set(gca,'XLim',[1 2],'XColor','none','YColor','none');
legend({"PV","SOM","EXC","PC 1","PC 2","PC 3"},'NumColumns',2,'box','off');



%% PLOT INSET FIGURE SHOWING TOP 3 PCs OF EACH CELL TYPE (DATASET 2)

targetRep = 1;  % target model iteration to visualize
pcVisNum = 3; % number of PCs to visualize
figure('theme','light','color',[1 1 1]);
tiledlayout(2,2);
nexttile(1)
spx = [-0.15 -0.05 0.05 0.15];
ntx = [1 2 4];

for nc = 1 : length(cellType)
    nexttile(ntx(nc)); hold on
    for np = 1 : pcVisNum
        locResps = PC.(cellType(nc)).score{targetRep}(:,np);
        meanYs = zeros(nFreq,1);
        for na = 1 : nAmp
            yVals = locResps(trainAmp == ampValAll(na));
            xVals = [ampValAll(na) ampValAll(na)]+spx(np);
            meanYs(na) = mean(yVals);
            plot(xVals,mkCI(yVals),'-','Color',colors.(cellType(nc))); % CI bar
            plot(xVals(1),meanYs(na),'.','Color',colors.(cellType(nc))) % mean point
        end 
        switch np
            case 1
                plot(ampValAll+spx(np),meanYs,'-','Color',colors.(cellType(nc)));
            case 2
                plot(ampValAll+spx(np),meanYs,'--','Color',colors.(cellType(nc)));
            case 3
                plot(ampValAll+spx(np),meanYs,':','Color',colors.(cellType(nc)));
        end
    end
    set(gca,'YLim',[-4 6],'YTick',-4:2:6,'XTick',ampValAll,'XLim',[-0.2 7.9],...
        'XTickLabelRotation',45,'FontName','Arial');
    ylabel("\DeltaF/F",'FontName','Arial');
    xlabel("Amplitude (\mum)",'FontName','Arial');
    hold off
end

% LEGEND IN SPOT 3
nexttile(3); hold on;
for nc = 1 : 3
plot(0,0,'.','Color',colors.(cellType(nc)),'MarkerSize',12);
end
plot([0 0],[0 0],'k-'); 
plot([0 0],[0 0],'k--');
plot([0 0],[0 0],'k:');
set(gca,'XLim',[1 2],'XColor','none','YColor','none');
legend({"PV","SOM","EXC","PC 1","PC 2","PC 3"},'NumColumns',2,'box','off',...
    'FontName','Arial','location','northeast');