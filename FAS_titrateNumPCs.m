% FAS_titrateNumPCs
%
% anotha one
%

cd('C:\Users\skich\Desktop\WORK');
load('modelAccuracy_dataset1_PCtiter.mat');

%% perform population decoding

nIter = 20; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
targetAmp = 5; % 0.33,0.76,1.6,3.5,7.7
usePCA = true; % use PCA step prior to model training\
useFreqEffNeu = true; % use neurons with significant anova freq effects
useAllNeu = false; % use all available neurons for modeling (no random cell draws)
cellType = ["PV","SOM","EXC"]; % order for plotting
nCellDraw = 280; % 280 dataset1, 160 datset2
nCellType = length(cellType);

% Create output structure for storing PCA information
PC = struct();
for nc = 1 : length(cellType)
    PC.(cellType(nc)) = struct();
end

nPCused = 1:1:15;

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
    trainFreq = repmat(freqVal,trainTrial,1);
    testFreq = repmat(freqVal,testTrial,1);
    predFreqAll = zeros(nIter,nCellType,nFreq*testTrial);
    testFreqAll = zeros(nIter,nCellType,nFreq*testTrial);
    LR_rsqs = zeros(nIter,nCellType);

    if useFreqEffNeu == true
        if multiAmp == true
            idx = find(Tc.identity == cellType(nc) & (is.FreqMod | is.DualMod | is.InterMod));
        elseif multiAmp == false
            idx = find(Tc.identity == cellType(nc) & Tc.selectivityPVal < pValThresh);
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
            if multiAmp
                cellAmp = Tc.amplitudes{selectedCells(nk)};
            end
                      
            for nl = 1:nFreq
                % randomly draw trials at each frequency, no replacement
                if multiAmp
                    trialNum = find(cellFreq == freqVal(nl) & cellAmp == ampValAll(targetAmp));
                else
                    trialNum = find(cellFreq == freqVal(nl));
                end

                selectedTrials = randsample(trialNum, trainTrial+testTrial, false);

                % Zscore data before splitting
                % selectedTrials = zscore(selectedTrials);

                % populate training and testing arrays
                trainResp(nl,:,nk) = cellResp(selectedTrials(1:trainTrial));
                testResp(nl,:,nk) = cellResp(selectedTrials(trainTrial+1:end));
 
            end
        end

        % Reshape matricies
        trainResp = reshape(trainResp, nFreq*trainTrial,nCellDraw);
        testResp = reshape(testResp, nFreq*testTrial,nCellDraw);

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
            mdl = fitcdiscr(trainResp_red,trainFreq,'DiscrimType','linear');
            predFreq = predict(mdl,testResp*coeff(:,nPCused(1):nPCused(np)));
    
            % collect output
            PCTall(ni,nc,np) = mean(predFreq == testFreq); %accuracy
            predFreqAll(ni,nc,:) = predFreq;
            testFreqAll(ni,nc,:) = testFreq;
            MDLall{ni,nc,np} = mdl;
        end
    end 

    % end %N PCS USED
    fprintf(newline)
end





%%  PLOT ACCURACY TOGETHER 
figure('Color',[1 1 1],'Theme','Light'); hold on;
% plot fake data for correct legend order
% plot([-1 -1],[-1 -1],'k-'); 
% plot([-1 -1],[-1 -1],'k--');
fill([-1 -1 -1], [-1 -1 -1],colors.PV,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.SOM,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.EXC,'FaceAlpha',0.5,'EdgeColor','none');
nPCused = 1:15;
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
        plot(nPCused,yMeans,'Marker','.','LineStyle','-','LineWidth',1,'Color',colors.(cellType(nc)));
    end
end

% plot([15 500],[20 20],'r--'); % chance performance level 
set(gca,'XLim',[nPCused(1)-0.02 nPCused(end)+0.5],'XTick',1:15,'YLim',[20 100],...
    'YTick',10:10:100);
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Number of PC predictors in model','FontName','Arial');
title('Frequency decoding performance by number of PCs','FontName','Arial');
subtitle('From 280 neurons with frequency selectivity','FontName','Arial');
% legend({"Freq. Mod.","All neurons","PV","SOM","EXC"},'Location','NorthWest','box','off');
legend({"PV","SOM","EXC"},'Location','West','box','off');


%%  PLOT ACCURACY TOGETHER (PC2-15)

% load("C:\Users\skich\Desktop\WORK\modelAccuracy_dataset1_PC2-15.mat")

figure('Color',[1 1 1],'Theme','Light'); hold on;
% plot fake data for correct legend order
% plot([-1 -1],[-1 -1],'k-'); 
% plot([-1 -1],[-1 -1],'k--');
fill([-1 -1 -1], [-1 -1 -1],colors.PV,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.SOM,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.EXC,'FaceAlpha',0.5,'EdgeColor','none');
nPCused = 2:15;


for nc = 1 : length(cellType)
    yMeans = mean(squeeze(PCTall(:,nc,:)),1).*100;
    yCIs = zeros(length(nPCused),2); % confidence intervals (95%)
    for np = 1 : length(nPCused)
        yCIs(np,:) = mkCI(squeeze(PCTall(:,nc,np))).*100;
    end
    xLooped = [nPCused flip(nPCused)]; 
    % yLooped = [yMeans+ySems flip(yMeans-ySems)]; when using SEM
    yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
    fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.5,'EdgeColor','none');
    plot(nPCused,yMeans,'Marker','.','LineStyle','-','LineWidth',1,'Color',colors.(cellType(nc)));
end


% plot([15 500],[20 20],'r--'); % chance performance level 
set(gca,'XLim',[nPCused(1)-0.02 nPCused(end)+0.5],'XTick',1:15,'YLim',[35 70],...
    'YTick',10:10:100);
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Number of PC predictors in model','FontName','Arial');
title('Frequency decoding performance by number of PCs','FontName','Arial');
subtitle('From 280 neurons with frequency selectivity','FontName','Arial');
% legend({"Freq. Mod.","All neurons","PV","SOM","EXC"},'Location','NorthWest','box','off');
legend({"PV","SOM","EXC"},'Location','West','box','off');


%%  PLOT ACCURACY TOGETHER (PC2-15) PV & EXC, (3-15) SOM AS DELTA VALUES

% cd('C:\Users\skich\Box\Tactile_Synchrony\2P-Data\Cleaned\datasets');
% load('modelAccuracy_perPC_dataset1_PC2-15_PC3-15.mat');

figure('Color',[1 1 1],'Theme','Light'); hold on;
% plot fake data for correct legend order
% plot([-1 -1],[-1 -1],'k-'); 
% plot([-1 -1],[-1 -1],'k--');
fill([-1 -1 -1], [-1 -1 -1],colors.PV,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.SOM,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.EXC,'FaceAlpha',0.5,'EdgeColor','none');
nPCused = 1:14;

nc = 3; % EXC NEURONS
nPCused = 1:14;
yMeans = mean(squeeze(PCTall(:,nc,:)),1).*100;
yMeanDelta = diff([0 yMeans]);
yCIs = zeros(length(nPCused),2); % confidence intervals (95%)
for np = 1 : length(nPCused)
    yCIs(np,:) = mkCI(squeeze(PCTall(:,nc,np))).*100;
end
yCIsZeroed =  (yCIs - repmat(yMeans',1,2)) + repmat(yMeanDelta',1,2);
xLooped = [nPCused flip(nPCused)]; 
% yLooped = [yMeans+ySems flip(yMeans-ySems)]; when using SEM
% yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
yLooped = [yCIsZeroed(:,1)' flip(yCIsZeroed(:,2))'];
fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.5,'EdgeColor','none');
plot(nPCused,yMeanDelta,'Marker','.','LineStyle','-','LineWidth',1,'Color',colors.(cellType(nc)));

nc = 2; % SOME NEURONS
nPCused = 1:13;
yMeans = mean(squeeze(PCTall3(:,nc,:)),1).*100;
yMeanDelta = diff([0 yMeans]);
yCIs = zeros(length(nPCused),2); % confidence intervals (95%)
for np = 1 : length(nPCused)
    yCIs(np,:) = mkCI(squeeze(PCTall3(:,nc,np))).*100;
end
yCIsZeroed =  (yCIs - repmat(yMeans',1,2)) + repmat(yMeanDelta',1,2);
xLooped = [nPCused flip(nPCused)]; 
% yLooped = [yMeans+ySems flip(yMeans-ySems)]; when using SEM
% yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
yLooped = [yCIsZeroed(:,1)' flip(yCIsZeroed(:,2))'];
fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.5,'EdgeColor','none');
plot(nPCused,yMeanDelta,'Marker','.','LineStyle','-','LineWidth',1,'Color',colors.(cellType(nc)));

nc = 1; % PV NEURONS
nPCused = 1:14;
yMeans = mean(squeeze(PCTall(:,nc,:)),1).*100;
yMeanDelta = diff([0 yMeans]);
yCIs = zeros(length(nPCused),2); % confidence intervals (95%)
for np = 1 : length(nPCused)
    yCIs(np,:) = mkCI(squeeze(PCTall(:,nc,np))).*100;
end
yCIsZeroed =  (yCIs - repmat(yMeans',1,2)) + repmat(yMeanDelta',1,2);
xLooped = [nPCused flip(nPCused)]; 
% yLooped = [yMeans+ySems flip(yMeans-ySems)]; when using SEM
% yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
yLooped = [yCIsZeroed(:,1)' flip(yCIsZeroed(:,2))'];
fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.5,'EdgeColor','none');
plot(nPCused,yMeanDelta,'Marker','.','LineStyle','-','LineWidth',1,'Color',colors.(cellType(nc)));


% plot([15 500],[20 20],'r--'); % chance performance level 
set(gca,'XLim',[0.5 14],'XTick',1:15,'YLim',[-5 40],...
    'YTick',0:10:100);
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Number of PC predictors in model','FontName','Arial');
title('Model performance gain from consecutive PC addition','FontName','Arial');
subtitle('PC2-15 (PV & EXC), PC3-15 (SOM), model size (n=270)','FontName','Arial');
% legend({"Freq. Mod.","All neurons","PV","SOM","EXC"},'Location','NorthWest','box','off');
legend({"PV","SOM","EXC"},'Location','West','box','off');



%%  PLOT ACCURACY TOGETHER AS DELTA VALUESsize(yCIs)

cd('C:\Users\skich\Desktop\WORK');
load("C:\Users\skich\Desktop\WORK\modelAccuracy_dataset1_PCtiter.mat")

nPCused = 1:15;
figure('Color',[1 1 1],'Theme','Light'); hold on;
% plot fake data for correct legend order
plot([-1 -1],[-1 -1],'k-'); 
% plot([-1 -1],[-1 -1],'k--');
fill([-1 -1 -1], [-1 -1 -1],colors.PV,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.SOM,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.EXC,'FaceAlpha',0.5,'EdgeColor','none');

% for nt = 1 : size(PCTall,4)
nt = 1;
for nc = 1 : length(cellType)
    yMeans = mean(squeeze(PCTall(:,nc,nPCused,nt)),1).*100;
    yDelta = [yMeans(1) diff(yMeans)];
    % ySems = sem(squeeze(PCTall(:,nc,:,nt)),1).*100; % standard error
    yCIs = zeros(length(nPCused),2); % confidence intervals (95%)
    for np = 1 : length(nPCused)
        yCIs(np,:) = mkCI(squeeze(PCTall(:,nc,np))).*100;
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
            plot(nPCused,yDelta,'Marker','.','MarkerSize',8,'LineStyle','-','LineWidth',1,'Color',...
                colors.(cellType(nc)));
        case 2
            % plot(nPCused,yMeans,'--','LineWidth',1,'Color',colors.(cellType(nc)));
            plot(nPCused,yDelta,'--','LineWidth',1,'Color',colors.(cellType(nc)));
    end
end
% end

% plot([15 500],[20 20],'r--'); % chance performance level 
% set(gca,'XLim',[1 40],'YLim',[60 100],'YTick',60:10:100);
% set(gca,'XLim',[1 6],'XTick',1:6,'YLim',[-5 80],'YTick',0:10:80);
% ylabel('model accuracy (%)','FontName','Arial');
ylabel('\Delta Frequency decoding accuracy (%)','FontName','Arial');
xlabel('PC number','FontName','Arial');
title('Model performance gain per PC','FontName','Arial');
subtitle('From 280 frequency selective neurons','FontName','Arial');
% legend({"Freq. Mod.","All neurons","PV","SOM","EXC"},'Location','NorthWest');
% legend({"Freq. Mod.","PV","SOM","EXC"},'Location','NorthWest','Box','off');
set(gca,'XLim',[0.9 15],'XTick',1:2:15,'YLim',[-5 80],'YTick',0:10:80,'Box','off');
% set(gca,'XLim',[2 20],'XTick',[2:2:20],'YLim',[-5 40],'YTick',0:10:80,'XScale','log','Box','off');


%% PLOT INSET FIGURE SHOWING TOP 3 PCs OF EACH CELL TYPE (DATASET 1)

targetRep = 2;  % target model iteration to visualize
pcVisNum = 3; % number of PCs to visualize
figure('theme','light','color',[1 1 1]);
tiledlayout(2,3,'TileSpacing','tight','Padding','tight');
nexttile(1)
spx = [-10 0 10];
ntx = [1 2 3];

for nc = 1 : length(cellType)
    nexttile(ntx(nc)); hold on
    if nc == 2
        mdl = MDLall{targetRep,nc,pcVisNum+1}; 
    else
        mdl = MDLall{targetRep,nc,pcVisNum}; 
    end
    PCresps = zeros(trainTrial,5,pcVisNum);
    for np = 1 : pcVisNum
        for nf = 1 : nFreq
            if nc == 2
                PCresps(:,nf,np) = mdl.X(mdl.Y == freqVal(nf),np+1);
            else
                PCresps(:,nf,np) = mdl.X(mdl.Y == freqVal(nf),np);
            end
            yMean = mean(PCresps(:,nf,np));
        end
    end

    for np = 1 : pcVisNum
        xVals = freqVal + spx(np); % frequency + xspacing  
        locResps = squeeze(PCresps(:,:,np));
        for nf = 1 : nFreq
            plot([xVals(nf) xVals(nf)],mkCI(locResps(:,nf)),'-','Color',colors.(cellType(nc))); % CI bar
            plot(xVals(nf),mean(locResps(:,nf)),'.','Color',colors.(cellType(nc))) % mean point
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
    % set(gca,'YLim',[-10 8],'YTick',-8:4:8,'XTick',freqVal,'XLim',[280 1120],'FontSize',15);
    set(gca,'YLim',[-5 5],'YTick',-4:2:4,'XTick',freqVal,'XLim',[280 1120],'FontSize',15);
    if nc == 1 
        ylabel("\DeltaF/F",'FontName','Arial');
    elseif nc == 2
        xlabel("Frequency (Hz)",'FontName','Arial');
    end
    hold off
end

% LEGEND IN SPOT 3
nexttile(4); hold on;
for nc = 1 : 3
plot(0,0,'.','Color',colors.(cellType(nc)),'MarkerSize',12);
end
plot([0 0],[0 0],'k-'); 
plot([0 0],[0 0],'k--');
plot([0 0],[0 0],'k:');
set(gca,'XLim',[1 2],'XColor','none','YColor','none','FontSize',12);
% legend({"PV","SOM","EXC","PC 1","PC 2","PC 3"},'NumColumns',2,'box','off','FontSize',15);
legend({"PV","SOM","EXC","PC 2","PC 3","PC 4"},'NumColumns',2,'box','off','FontSize',15);


%% PLOT SHOWING PC CORRELATIONS WITH STIM AMPLITUDE COMPONENT
% cd('C:\Users\skich\Desktop\WORK');
% load("C:\Users\skich\Desktop\WORK\modelResults_dataset1_PC1-15.mat");

% MDLall
% Reps (200) x CellType (3) x nPCs (15 or 14)
ampValsOg = [20; 35; 70; 45; 30];
targetBin = 15;  % target bin
numPC = length(MDLall{end}.PredictorNames);
numReps = size(MDLall,1);
allCorrs = zeros(numReps,length(cellType),numPC); % Reps x CellType x PCnum
xSP = [-0.2 0 0.2];

for nc = 1 : length(cellType) % loop cell type
    for nr = 1 : numReps % loop model reps
        mdl = MDLall{nr,nc,targetBin}; 
        numTrials = length(mdl.Y) / length(freqVal);
        PCresps = zeros(numTrials,length(freqVal),numPC);
        for np = 1 : numPC
            for nf = 1 : nFreq
                PCresps(:,nf,np) = mdl.X(mdl.Y == freqVal(nf),np);
                yMean = mean(PCresps(:,nf,np));
            end
        end
        locMeanResps = squeeze(mean(PCresps,1));
        for np = 1 : numPC
            allCorrs(nr,nc,np) = corr(locMeanResps(:,np),ampValsOg,'Type','Spearman');
        end
    end
end

% NEW WAY (2025-12-17)
figure('theme','light','color',[1 1 1]); hold on;
for nc = 1 : length(cellType)
    allCIs = zeros(numPC,2);
    for np = 1 : numPC
        allCIs(np,:) = mkCI(allCorrs(:,nc,np));
    end
    % SHADED ERROR PLOT
    xVals = [1:numPC flip(1:numPC)]';
    yVals = cat(1,allCIs(:,1),flip(allCIs(:,2)));
    fill(xVals,yVals,colors.(cellType(nc)),'FaceAlpha',0.5,'EdgeColor','none');
    % PLOT LINE BETWEEN POINTS 
    plot(1:numPC,squeeze(mean(allCorrs(:,nc,:))),'-','Color',colors.(cellType(nc)));
    plot(1:numPC,squeeze(mean(allCorrs(:,nc,:))),'.','MarkerSize',10,'Color',colors.(cellType(nc)));
end

% OLD WAY (2025-12-16)
% figure('theme','light','color',[1 1 1]); hold on;
% for nc = 1 : length(cellType)
%     % PLOT LINE BETWEEN POINTS 
%     plot((1:numPC)+xSP(nc),squeeze(mean(allCorrs(:,nc,:))),'-','Color',colors.(cellType(nc)))
% end
% 
% for np = 1 : numPC
%     for nc = 1 : length(cellType)
%         % plot(repmat(np+xSP(nc),2,1),mkCI(allCorrs(:,nc,np)),'-','LineWidth',2,'Color',colors.(cellType(nc)));
%         % plot(np+xSP(nc),locMean,'k.','MarkerSize',12);
%         locMean = mean(allCorrs(:,nc,np));
%         locCIs = mkCI(allCorrs(:,nc,np)) - locMean;
%         errorbar(np+xSP(nc),locMean,locCIs(1),locCIs(2),'LineWidth',1.2,'Color',colors.(cellType(nc)));
%         plot(np+xSP(nc),locMean,'.','MarkerSize',10,'Color',colors.(cellType(nc)));
%     end
% end
set(gca,'YLim',[-0.6 1.01],'XTick',[1 5 10 15],'FontName','Arial');
xlabel('PC Number','FontName','Arial');
ylabel('Avg Correlation (Spearman)','FontName','Arial');
title('PC Correlation with stimulus amplitude','FontName','Arial');
yline(0,':');
legend({'','PV','','','SOM','','','EXC','',''},"Box","off",'FontName','Arial','Location','Southeast');







%% PLOT INSET FIGURE SHOWING TOP 3 PCs OF EACH CELL TYPE (DATASET 2)

targetRep = 1;  % target model iteration to visualize
pcVisNum = 4; % number of PCs to visualize
figure('theme','light','color',[1 1 1]);
tiledlayout(2,2);
nexttile(1)
spx = [-15 -5 5 15];
ntx = [1 2 4];

for nc = 1 : length(cellType)
    nexttile(ntx(nc)); hold on
    for np = 1 : pcVisNum
        locResps = PC.(cellType(nc)).score{targetRep}(:,np);
        meanYs = zeros(nFreq,1);
        for nf = 1 : nFreq
            yVals = locResps(trainFreq == freqVal(nf));
            xVals = [freqVal(nf) freqVal(nf)]+spx(np);
            meanYs(nf) = mean(yVals);
            plot(xVals,mkCI(yVals),'-','Color',colors.(cellType(nc))); % CI bar
            plot(xVals(1),meanYs(nf),'.','Color',colors.(cellType(nc))) % mean point
        end 
        switch np
            case 1
                plot(freqVal+spx(np),meanYs,'-','Color',colors.(cellType(nc)));
            case 2
                plot(freqVal+spx(np),meanYs,'--','Color',colors.(cellType(nc)));
            case 3
                plot(freqVal+spx(np),meanYs,':','Color',colors.(cellType(nc)));
            case 4
                plot(freqVal+spx(np),meanYs,'-.','Color',colors.(cellType(nc)));
        end
    end
    set(gca,'YLim',[-10 8],'YTick',-8:2:8,'XTick',freqVal,'XLim',[280 1620]);
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
plot([0 0],[0 0],'k-.');
set(gca,'XLim',[1 2],'XColor','none','YColor','none');
legend({"PV","SOM","EXC","PC 1","PC 2","PC 3","PC 4"},'NumColumns',2);



%% MEASURE VARIANCE EXPLAINED FOR EACH CELL TYPE

targetRep = 3;  % target model iteration to visualize
figure('theme','light','color',[1 1 1]); hold on;

for nc = 1 : length(cellType)
    latent = PC.(cellType(nc)).latent{targetRep};
    totalVariance = sum(latent);
    CVA = cumsum((latent./totalVariance).*100); %cumulative variance explained
    plot(CVA,'LineStyle','-','Color',colors.(cellType(nc)),'Marker','.');
end

ylabel('% Variance Explained','FontName','Arial');
xlabel('PC number','FontName','Arial');
title('Cumulative variance explained per cell type ','FontName','Arial');
set(gca,'XLim',[0 100],'YLim',[0 100],'FontName','Arial');