% FAS_titrateSelectivity
%
%

% Loop through each cell type
% --> Loop through iterations
% ----> Draw 280 random neurons from freq. sel. pool
% ----> Sort 280 neurons by selectivity strength
% ----> Loop through sel bin to draw from 50:10:280



%% perform population decoding

nIter = 200; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
usePCA = true; % use PCA step prior to model training
kStart = 1; %first PC to include in training / testing
kComp = 15; %number of PCs to use in model training / testing
useFreqEffNeu = true; % use neurons with significant anova freq effects
useAllNeu = false; % use all available neurons for modeling (no random cell draws)
cellType = ["PV","SOM","EXC"]; % order for plotting
nCellDraw = 280;
binSize = 50; % number of neurons per selectivity bin
binStepSize = 10; %number of neurons to subtract / add between each bin


% Create output structure for storing PCA information
PC = struct();
for nc = 1 : length(cellType)
    PC.(cellType(nc)) = struct();
end

if multiAmp == true
    nSelBins = 20:20:1000;
elseif multiAmp == false
    startIndx = 1:binStepSize:(nCellDraw-binSize+1); % 1:10:231
    stopIndx = binSize:binStepSize:nCellDraw; % 50:10:280
    nSelBins = cat(2,startIndx',stopIndx');
end

% CREATE OUTPUT ARRAYS
% Performance of Cell Types (all)
% PCTall = []; %for populating with accuracy measurements
PCTall = zeros(nIter,nCellType,length(nSelBins)); % (nj,ni,nd): 200x3x24
% Models (all)
MDLall = cell(nIter,nCellType,length(nSelBins)); % (nj,ni,nd)

for nc = 1:length(cellType)% 
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
       
        % calculate selectivity metric for cells in selectedCells
        selectedCells(:,2) = 0; % store selectivity values in second colum
        for ns = 1 : length(selectedCells)
            locResps = Tc.responses{selectedCells(ns)};
            locFreqs = Tc.frequencies{selectedCells(ns)};
            meanResps = zeros(nFreq,1);
            for nf = 1 : nFreq
                meanResps(nf) = mean(locResps(locFreqs == freqVal(nf)));
            end
            selectedCells(ns,2) = quickSelectivity(meanResps);
        end
        [~,idx1] = sort(selectedCells(:,2),'descend'); % sort by descending selectivity strength
        strongestCells = selectedCells(idx1,1); % store sorted cell numbers in sorted order


        for nb = 1 : length(nSelBins) % 24 

            binIndx = nSelBins(nb,1):nSelBins(nb,2);

            % pre-allocate training and testing data array
            trainResp = zeros(nFreq,trainTrial,binSize);
            testResp = zeros(nFreq,testTrial,binSize);
    
            % loop over selected cells
            for nk = 1:binSize %1:50

                cellResp = Tc.responses{strongestCells(binIndx(nk))};
                cellFreq = Tc.frequencies{strongestCells(binIndx(nk))};
                          
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
            trainResp = reshape(trainResp, nFreq*trainTrial,binSize);
            testResp = reshape(testResp, nFreq*testTrial,binSize);
    
            % Z-Scoring 
            for nz = 1 : binSize
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
    
    
            % kComp = floor(trainTrial*nFreq*0.5);
            trainResp_red = score(:,kStart:kComp);
            mdl = fitcdiscr(trainResp_red,trainFreq,'DiscrimType','linear');
            predFreq = predict(mdl,testResp*coeff(:,kStart:kComp));
    
            % collect output
            PCTall(ni,nc,nb) = mean(predFreq == testFreq); %accuracy
            predFreqAll(ni,nc,:) = predFreq;
            testFreqAll(ni,nc,:) = testFreq;
            MDLall{ni,nc,nb} = mdl;

        end

    end 

    % end %N PCS USED
    fprintf(newline)
end




%%  PLOT ACCURACY TOGETHER 

figure('Color',[1 1 1],'Theme','Light'); hold on;
% plot fake data for correct legend order
fill([-1 -1 -1], [-1 -1 -1],colors.PV,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.SOM,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.EXC,'FaceAlpha',0.5,'EdgeColor','none');

for nc = 1 : length(cellType)
    yMeans = mean(squeeze(PCTall(:,nc,:)),1).*100;
    % ySems = sem(squeeze(PCTall(:,nc,:,nt)),1).*100; % standard error
    yCIs = zeros(length(nSelBins),2); % confidence intervals (95%)
    for nb = 1 : length(nSelBins)
        yCIs(nb,:) = mkCI(squeeze(PCTall(:,nc,nb))).*100;
    end
    xVals = 1:length(nSelBins);
    xLooped = [xVals flip(xVals)]; 
    % yLooped = [yMeans+ySems flip(yMeans-ySems)]; when using SEM
    yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
    fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.5,'EdgeColor','none');
    plot(xVals,yMeans,'-','LineWidth',1,'Color',colors.(cellType(nc)));
end


% plot([15 500],[20 20],'r--'); % chance performance level 
set(gca,'XLim',[0 25],'YLim',[45 75],'YTick',10:10:100);
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Selectivity Bin (descending strength)','FontName','Arial');
title('Frequency decoding performance by selectivity strength','FontName','Arial');
subtitle('From 50 neuron bins','FontName','Arial');
legend({"PV","SOM","EXC"},'Location','NorthWest');

























