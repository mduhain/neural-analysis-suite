% FAS_indvMouseOmit.m
%
% sub script from focusedAnalysisScript.m
%
% perform training by omiting individual mice, asess deviance from baseline
% accuracy

%% perform population decoding

% TUNABLE INPUTS
nIter = 400;   % Number of model iterations to run
trainPortion = 0.8;   % Percentage of trials to be used in training (remaining for testing)
targetAmp = 5;   % use amplitude: [0.33, 0.76, 1.6, 3.5, 7.7]
usePCA = true;   % use PCA step prior to model training
kComp = 15;   %number of PCs to use in model training / testing
useFreqEffNeu = true;   % use neurons with significant anova freq effects
useAllNeu = false;   % use all available neurons for modeling
includeErrBar = true;   % include error bars in plots
cellType = ["PV"; "SOM"; "EXC"];   % order for plotting
nCellDraws = 200;   % Number of neurons to draw from each pool

% PRE-POPULATE ARRAYS
PCTall = [];   %for populating with accuracy measurements
PerfCellType = zeros(nIter,nCellType,length(mouseNums)+1); % performance array

for nd = 1 : length(mouseNums) + 1

    if nd == 1 
        disp("All mice: ");
    else
        fprintf(newline);
        disp(strcat("no ",mouseNums(nd-1)));
    end

    % PRE-ALLOCATE ARRAYS
    trainTrial = round(minTrialNum*trainPortion);
    testTrial = round(minTrialNum*(1-trainPortion));
    trainFreq = repmat(freqVal,trainTrial,1);
    testFreq = repmat(freqVal,testTrial,1);
    predFreqAll = zeros(nIter,nCellType,nFreq*testTrial);
    testFreqAll = zeros(nIter,nCellType,nFreq*testTrial);
    LR_rsqs = zeros(nIter,nCellType);

    for ni = 1:length(cellType)
        fprintf(strcat(" ",cellType(ni)))
        
        % select cell type
        if useFreqEffNeu == true
            if multiAmp == true
                idx = find(Tc.identity == cellType(ni) & is.InterMod);
            elseif multiAmp == false
                if nd == 1 % First rep, no cell removal
                    idx = find(Tc.identity == cellType(ni) & is.selective);
                else
                    idx = find(Tc.identity == cellType(ni) & is.selective & ~is.mouseNum(:,nd-1));
                end
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
                    if multiAmp == true
                        trialNum = find(cellFreq == freqVal(nl) & cellAmp == ampValAll(targetAmp));
                    elseif multiAmp == false
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
    
            % turn array into matrices and zscore each cell
            trainResp = reshape(trainResp, nFreq*trainTrial,nCellDraw);
            %trainResp = normalize(trainResp,1,"zscore");
    
            testResp = reshape(testResp, nFreq*testTrial,nCellDraw);
            %testResp = normalize(testResp,1,"zscore");
    
            % alternate zscoring approach
            for nz = 1 : nCellDraw
                mu = mean(trainResp(:,nz));
                sigma = std(trainResp(:,nz));
                % Z-score training
                trainResp(:,nz) = (trainResp(:,nz) - mu) ./ sigma;
                % Z-score test using training stats
                testResp(:,nz) = (testResp(:,nz) - mu) ./ sigma;
            end
    
            if usePCA
                % WITH PCA 
                [coeff,score] = pca(trainResp);
                % kComp = floor(trainTrial*nFreq*0.5);
                trainResp_red = score(:,1:kComp);
                mdl = fitcdiscr(trainResp_red,trainFreq,'DiscrimType','linear');
                predFreq = predict(mdl,testResp*coeff(:,1:kComp));
            else
                % WITHOUT PCA
                mdl = fitcdiscr(trainResp,trainFreq,'DiscrimType','linear');
                predFreq = predict(mdl,testResp);
            end

            % collect output
            acc = mean(predFreq == testFreq);
            predFreqAll(nj,ni,:) = predFreq;
            testFreqAll(nj,ni,:) = testFreq;
            PerfCellType(nj,ni,nd) = acc;
    
        end
    end
    fprintf(newline);
end

% figure('Theme','light'); hold on;
% for nd = 1 : length(mouseNums)+1
%     xVals = (1:3)+3*(nd-1);
%     yVals = squeeze(PerfCellType(:,:,nd));
%     v = violinplot(xVals,yVals);
%     v(1).FaceColor = colors.(cellType(1));
%     v(2).FaceColor = colors.(cellType(2));
%     v(3).FaceColor = colors.(cellType(3));
% end

% Change in model performance chart
figure('Theme','light','Color',[1 1 1],'Position',[100 100 500 500]);
tiledlayout(3,1)
for nc = 1 : length(cellType)
    yVals = squeeze(PerfCellType(:,nc,:));
    yMeans = mean(yVals);
    yBaselineAvg = yMeans(1);
    yDelta = yMeans(2:end) - yBaselineAvg;
    ySems = arrayfun(@(col) sem(yVals(:,col)), 1:size(yVals,2));
    nexttile(nc); hold on;
    plot([0 length(yDelta)+1],[0 0],'k--');
    errorbar(yDelta,ySems(2:end),"LineStyle","none",'Color',colors.(cellType(nc)))
    plot(yDelta,'.','Color',colors.(cellType(nc)),'MarkerSize',10);
    % STATS
    for nr = 1 : length(yDelta)
        [p,h] = ranksum(yVals(:,1),yVals(:,nr+1));
        if p < pValThresh
            plot(nr+0.2,yDelta(nr),'r*');
        end
    end
    % LABELS
    switch nc
        case 1 
            title("Individual mouse effects on model performance",'FontName','Arial');
            set(gca,'XColor','none','XLim',[0 16]);
        case 2
            ylabel('\Delta Performance','FontName','Arial');
            set(gca,'XColor','none','XLim',[0 16]);
        case 3
            xlabel('Mouse Number','FontName','Arial');
            set(gca,'XLim',[0 16]);
    end
    hold off;
end










