% FAS_selTypePerf_mixCellTypes.m
%
% Selectivity Type Performance from all cells mixed together
%
% Dataset2, built from FAS_selTypePerf_allCells.m
%
%



%% perform baseline population decoding (one per cell type)

nIter = 500; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
targetAmp = 5; % 0.33,0.76,1.6,3.5,7.7
kStart = 1; % first PC to include
kComp = 15; % total number of PCs to use in model training / testing
cellType = ["PV","SOM","EXC"]; % order for plotting
selType = "InterMod"; 
nCellDraws = 40; % over-estimation to poulate arrays

% CREATE OUTPUT ARRAYS
PCTall1 = zeros(nIter,length(cellType)); % (nj,nc)
MDLall1 = cell(nIter,length(cellType)); % (nj,nc)

for nc = 1:length(cellType) %nt = 1 : length(selType)
    disp(cellType(nc));
    trainTrial = round(minTrialNum*trainPortion);
    testTrial = round(minTrialNum*(1-trainPortion));
    trainFreq = repmat(freqVal,trainTrial,1);
    testFreq = repmat(freqVal,testTrial,1);
    tic;

    % find cells of this type
    % idx = find(Tc.identity == cellType(nc) & is.(selType));
    % idx = find(Tc.identity == cellType(nc)); % any neurons
    idx = find(Tc.identity == cellType(nc) & (is.FreqMod | is.DualMod | is.InterMod)); % any frequency effected neurons

    % reiterate random draws
    for nj = 1:nIter

        % randomly draw cells without replacement
        selectedCells = randsample(idx, nCellDraws, false);
       
        % pre-allocate training and testing data array
        trainResp = zeros(nFreq,trainTrial,nCellDraws);
        testResp = zeros(nFreq,testTrial,nCellDraws);

        % loop over selected cells
        for nk = 1:length(selectedCells)
            cellResp = Tc.responses{selectedCells(nk)};
            cellFreq = Tc.frequencies{selectedCells(nk)};
            cellAmp = Tc.amplitudes{selectedCells(nk)};
                      
            for nl = 1:nFreq
                % randomly draw trials at each frequency, no replacement
                trialNum = find(cellFreq == freqVal(nl) & cellAmp == ampValAll(targetAmp));
                selectedTrials = randsample(trialNum, trainTrial+testTrial, false);

                % populate training and testing arrays
                trainResp(nl,:,nk) = cellResp(selectedTrials(1:trainTrial));
                testResp(nl,:,nk) = cellResp(selectedTrials(trainTrial+1:end));
            end
        end

        % Reshape matricies
        trainResp = reshape(trainResp, nFreq*trainTrial,nCellDraws);
        testResp = reshape(testResp, nFreq*testTrial,nCellDraws);

        % Z-Scoring 
        for nz = 1 : nCellDraws
            mu = mean(trainResp(:,nz));
            sigma = std(trainResp(:,nz));
            trainResp(:,nz) = (trainResp(:,nz) - mu) ./ sigma;
            testResp(:,nz) = (testResp(:,nz) - mu) ./ sigma;
        end

        % Training
        [coeff,score,latent] = pca(trainResp);
        trainResp_red = score(:,kStart:kComp);
        mdl = fitcdiscr(trainResp_red,trainFreq,'DiscrimType','linear'); %train model
       
        predFreq = predict(mdl,testResp*coeff(:,kStart:kComp));

        % collect output
        PCTall1(nj,nc) = mean(predFreq == testFreq); % accuracy
        MDLall1{nj,nc} = mdl; % model parameters
    end
end



% PERFORM MIXED POPULATION DECODING
cellTypeX = ["EXC_PV","EXC_SOM"]; % order for plotting
cellRatio = [0.5 0.5; 0.7 0.3; 0.9 0.1];

% CREATE OUTPUT ARRAYS
PCTall2 = zeros(nIter,length(cellRatio),length(cellTypeX)); % (nj,nr,nc)
MDLall2 = cell(nIter,length(cellRatio),length(cellTypeX)); % (nj,nr,nc)

for nc = 1:length(cellTypeX) %nt = 1 : length(selType)
    disp(cellTypeX(nc));
    trainTrial = round(minTrialNum*trainPortion);
    testTrial = round(minTrialNum*(1-trainPortion));
    trainFreq = repmat(freqVal,trainTrial,1);
    testFreq = repmat(freqVal,testTrial,1);
    tic;

    % find cells of this type
    cellTargs = strsplit(cellTypeX(nc),'_');
    % idx1 = find(Tc.identity == cellTargs(1) & is.(selType));
    % idx2 = find(Tc.identity == cellTargs(2) & is.(selType));

    idx1 = find(Tc.identity == cellTargs(1) & (is.FreqMod | is.DualMod | is.InterMod));
    idx2 = find(Tc.identity == cellTargs(2) & (is.FreqMod | is.DualMod | is.InterMod));

    % loop through EXC:INH ratios
    for nr = 1 : length(cellRatio)
        nCellDraw1 = nCellDraws * cellRatio(nr,1);
        nCellDraw2 = nCellDraws * cellRatio(nr,2);

        % reiterate random draws
        for nj = 1 : nIter
    
            % randomly draw cells without replacement combine each group
            selectedCells = cat(1,randsample(idx1,nCellDraw1,false), randsample(idx2,nCellDraw2,false));
           
            % pre-allocate training and testing data array
            trainResp = zeros(nFreq,trainTrial,nCellDraws);
            testResp = zeros(nFreq,testTrial,nCellDraws);
    
            % loop over selected cells
            for nk = 1:length(selectedCells)
                cellResp = Tc.responses{selectedCells(nk)};
                cellFreq = Tc.frequencies{selectedCells(nk)};
                cellAmp = Tc.amplitudes{selectedCells(nk)};
                          
                for nl = 1:nFreq
                    % randomly draw trials at each frequency, no replacement
                    trialNum = find(cellFreq == freqVal(nl) & cellAmp == ampValAll(targetAmp));
                    selectedTrials = randsample(trialNum, trainTrial+testTrial, false);
    
                    % populate training and testing arrays
                    trainResp(nl,:,nk) = cellResp(selectedTrials(1:trainTrial));
                    testResp(nl,:,nk) = cellResp(selectedTrials(trainTrial+1:end));
                end
            end
    
            % Reshape matricies
            trainResp = reshape(trainResp, nFreq*trainTrial,nCellDraws);
            testResp = reshape(testResp, nFreq*testTrial,nCellDraws);
    
            % Z-Scoring 
            for nz = 1 : nCellDraws
                mu = mean(trainResp(:,nz));
                sigma = std(trainResp(:,nz));
                trainResp(:,nz) = (trainResp(:,nz) - mu) ./ sigma;
                testResp(:,nz) = (testResp(:,nz) - mu) ./ sigma;
            end
    
            % Training
            [coeff,score,latent] = pca(trainResp);
            trainResp_red = score(:,kStart:kComp);
            mdl = fitcdiscr(trainResp_red,trainFreq,'DiscrimType','linear'); %train model
           
            predFreq = predict(mdl,testResp*coeff(:,kStart:kComp));
    
            % collect output
            PCTall2(nj,nr,nc) = mean(predFreq == testFreq); % accuracy
            MDLall2{nj,nr,nc} = mdl; % model parameters
        end
    end
end




% PLOTTING
% ACCURACY AT N=40, EACH OF 3 CELL TYPES ALONE + COMBINATIONS

% blended colors
colors.EXC_PV = [0 0 0.6];
colors.EXC_SOM = [0.543 0.269 0.074];

figure('theme','light','color',[1 1 1]); hold on;
plot(-1,-1,'.','MarkerSize',8,'Color',colors.PV);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.SOM);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.EXC);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.EXC_PV);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.EXC_SOM);

% BASELINE PERFORMANCE
for nc = 1 : length(cellType)
    yVals = PCTall1(:,nc);
    yMean = mean(yVals)*100; % mean y value converted to percentage
    yCIs = (mkCI(yVals) - mean(yVals)).*100; % 95% confidence intervals centered at mean
    plot(nc,yMean,'.','MarkerSize',8,'Color',colors.(cellType(nc)));
    errorbar(nc,yMean,yCIs(1),yCIs(2),'Color',colors.(cellType(nc)));
end

% BASELINE PERFORMANCE
xVal = 4;
for nc = 1 : length(cellTypeX)
    for nr = 1 : length(cellRatio)
        yVals = PCTall2(:,nr,nc);
        yMean = mean(yVals)*100; % mean y value converted to percentage
        yCIs = (mkCI(yVals) - mean(yVals)).*100; % 95% confidence intervals centered at mean
        plot(xVal,yMean,'.','MarkerSize',8,'Color',colors.(cellTypeX(nc)));
        errorbar(xVal,yMean,yCIs(1),yCIs(2),'Color',colors.(cellTypeX(nc)));
        
        ratioTxt = [num2str(cellRatio(nr,1)*100) ':' num2str(cellRatio(nr,2)*100)];
        text(xVal,yMean+yCIs(2)+0.5,ratioTxt,"VerticalAlignment","middle",'Rotation',45,...
            'FontName','Arial','FontSize',10);
        xVal = xVal + 1;
    end
end

ylabel('model accuracy (%)','FontName','Arial');
xlabel('Cell type category','FontName','Arial');
title('Frequency decoding performance by combination','FontName','Arial');
subtitle(strcat("Amplitude: ",num2str(ampValAll(targetAmp),2),"\mum, model sizes: ",...
    num2str(nCellDraws)," neurons ea."),'FontName','Arial');
set(gca,'XLim',[0.5 9.5],'XTick',[],'FontName','Arial');
legend({"PV","SOM","EXC","EXC+PV","EXC+SOM"},'Location','Northeast','FontName','Arial','box','off');
