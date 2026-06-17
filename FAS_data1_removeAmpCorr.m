% FAS_data1_removeAmpCorr
%
%


%% perform population decoding

nIter = 200; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
targetAmp = 5; % 0.33,0.76,1.6,3.5,7.7
usePCA = true; % use PCA step prior to model training
kStart = 1; % first PC to include
kComp = 15; % number of PCs to use in model training / testing
useFreqEffNeu = true; % use neurons with significant anova freq effects
useAllNeu = false; % use all available neurons for modeling (no random cell draws)
nCellDraw = 150; %min available, across all corrThresh, across all cell types

corrThresh = flip(0.8:0.1:1);

% CREATE OUTPUT ARRAYS
% Performance of Cell Types (all)
% PCTall = []; %for populating with accuracy measurements
PCTall = zeros(nIter,nCellType,length(corrThresh)); % (nj,ni,nd)
% Models (all)
MDLall = cell(nIter,nCellType,length(corrThresh)); % (nj,ni,nd)

for nd = 1 : length(corrThresh)
    disp(strcat("val: ",num2str(corrThresh(nd))))
    trainTrial = round(minTrialNum*trainPortion);
    testTrial = round(minTrialNum*(1-trainPortion));
    trainFreq = repmat(freqVal,trainTrial,1);
    testFreq = repmat(freqVal,testTrial,1);
    predFreqAll = zeros(nIter,nCellType,nFreq*testTrial);
    testFreqAll = zeros(nIter,nCellType,nFreq*testTrial);
    LR_rsqs = zeros(nIter,nCellType);
    
    tic;
    for ni = 1:length(cellType)

        fprintf(strcat(" ",cellType(ni)))
    
        % select cell type
        if useFreqEffNeu == true
            if multiAmp == true
                idx = find(Tc.identity == cellType(ni) & (is.FreqMod | is.DualMod | is.InterMod));
            elseif multiAmp == false
                idx = find(is.(cellType(ni)) & is.selective & corrValues<corrThresh(nd));
                fprintf(strcat(" ",num2str(length(idx))));
            end
        elseif useFreqEffNeu == false
            idx = find(Tc.identity == cellType(ni));
        end
        
        % reiterate random draws
        for nj = 1:nIter
            
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
    
            % Training
            if usePCA
                % WITH PCA 
                [coeff,score] = pca(trainResp);
                % kComp = floor(trainTrial*nFreq*0.5);
                trainResp_red = score(:,kStart:kComp);
                mdl = fitcdiscr(trainResp_red,trainFreq,'DiscrimType','linear');
                predFreq = predict(mdl,testResp*coeff(:,kStart:kComp));
            else
                % WITHOUT PCA
                mdl = fitcdiscr(trainResp,trainFreq,'DiscrimType','linear');
                predFreq = predict(mdl,testResp);
            end
    
            % collect output
            PCTall(nj,ni,nd) = mean(predFreq == testFreq); %accuracy
            predFreqAll(nj,ni,:) = predFreq;
            testFreqAll(nj,ni,:) = testFreq;
            MDLall{nj,ni,nd} = mdl;
        end 
    end
    fprintf(newline)
end


%%  PLOT ACCURACY TOGETHER (NEW WAY?)

nCellDraws = 1:21;
figure('Color',[1 1 1],'Theme','Light'); hold on;

for nc = 1 : length(cellType)
    locMat = squeeze(PCTall(:,nc,:));
    yMeans = mean(locMat,1).*100;
    yCIs = zeros(size(locMat,2),2); % confidence intervals (95%)
    for nd = 1 : size(locMat,2)
        yCIs(nd,:) = mkCI(squeeze(PCTall(:,nc,nd))).*100;
    end
    xLooped = [nCellDraws flip(nCellDraws)]; 
    yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
    fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.5,'EdgeColor','none');
    plot(nCellDraws,yMeans,'-','LineWidth',1,'Color',colors.(cellType(nc)));
end

% plot([15 500],[20 20],'r--'); % chance performance level 
set(gca,'XLim',[1 21],'XTick',[1 5 10 15 20],'XTickLabel',...
    {'all','<0.95','<0.9','<0.85','<0.8'},'YLim',[75 90],'YTick',80:5:90);
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Neurons used','FontName','Arial');
title('Frequency decoding performance by cell type','FontName','Arial');
subtitle('As amplitude correlated neurons are removed','FontName','Arial');
legend({"PV","","SOM","","EXC"},'Location','NorthEast')