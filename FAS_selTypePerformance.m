% FOCUSED ANALYSIS SCRIPT
%
% Selectivity Type Performance
%
% Dataset1, selective draws vs. random draws
%
%

%% Preload model results

cd('C:\Users\skich\Desktop');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [5 amps, 200 iterations, 4 sel types, 15:5:N draws, 3 cell types]
% tic; load("C:\Users\skich\Desktop\modelResults_PCTall_MDLall.mat"); toc;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [5 amps, 200 iterations, 1 sel type (all freq), 15:5:N draws, 3 cell types]
% tic; load("C:\Users\skich\Desktop\modelResults_PCTall2_MDLall2.mat"); toc;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [top amp, 200 iterations, 4 sel types, 15:5:990 draws, 3EXC]
% tic; load("C:\Users\skich\Desktop\modelResults_PCTall_MDLall.mat"); toc;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%% perform population decoding

nIter = 200; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
usePCA = true; % use PCA step prior to model training
kStart = 1; % first PC to include
kComp = 15; % total number of PCs to use in model training / testing
useAllNeu = false; % use all available neurons per cell type for modeling
cellType = ["PV","SOM","EXC"]; % order for plotting
% cellType = "SOM"; % order for plotting
is.any = true(length(is.selective),1);
selType = ["selective","any"]; 
% selType = ["selective"]; 
neuStepSize = 15; % num neurons to add to model each step
nCellDraws = 33; % over-estimation to poulate arrays
maxDraw = 500;

% CREATE OUTPUT ARRAYS
PCTall = zeros(nIter,length(selType),nCellDraws,length(cellType)); % (nj,nt,nd,nc)
MDLall = cell(nIter,length(selType),nCellDraws,length(cellType)); % (nj,nt,nd,nc)
CVEall = zeros(nIter,length(selType),nCellDraws,length(cellType));
PCFall = cell(nIter,length(selType),nCellDraws,length(cellType)); % (nj,nt,nd,nc)

% % Create output structure for storing PCA information
% PC = struct();
% for nc = 1 : length(cellType)
%     PC.(cellType(nc)) = struct();
% end


for nc = 1:length(cellType) %nt = 1 : length(selType)
disp(cellType(nc));
% disp(selType(nt));

% IGNORE SOME PCs WITH SPECIFIC NEURONS
if nc == 2
    kStart = 1;
else
    kStart = 1;
end

for nt = 1 : length(selType) %nd = 1 : length(nCellDraws)
    fprintf(strcat(selType(nt)," "))
    %disp(strcat("Draw ",num2str(nCellDraws(nd))))

    if strcmp(cellType(nc),'EXC')
        maxNeuAvail = maxDraw;
    else % PV | SOM
        maxNeuAvail = sum(is.(cellType(nc)) & is.(selType(nt)));
        if maxNeuAvail > maxDraw
            maxNeuAvail = maxDraw;
        end
    end

    nCellDraws = kComp : neuStepSize : maxNeuAvail;
    trainTrial = round(minTrialNum*trainPortion);
    testTrial = round(minTrialNum*(1-trainPortion));
    trainFreq = repmat(freqVal,trainTrial,1);
    testFreq = repmat(freqVal,testTrial,1);
    
    for nd = 1 : length(nCellDraws) %nc = 1:length(cellType)
        %fprintf(strcat(" ",cellType(nc)))
        fprintf(strcat(num2str(nCellDraws(nd))," "));
        nCellDraw = nCellDraws(nd);
        idx = find(Tc.identity == cellType(nc) & is.(selType(nt)));
        
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
            [coeff,score,latent] = pca(trainResp);
            % PC.(cellType(nc)).coeff{ni} = coeff;
            % PC.(cellType(nc)).score{ni} = score;
            % PC.(cellType(nc)).latent{ni} = latent;

            % Percent variance explained
            explained = 100 * latent / sum(latent);
            CVEall(nj,nt,nd,nc) = sum(explained(1:15));
            
            % % Bar plot + cumulative line
            % figure;
            % bar(explained,'FaceColor',[0.6 0.6 0.9]); hold on;
            % plot(cumsum(explained),'r-o','LineWidth',1.5);
            % xlabel('Principal Component');
            % ylabel('Variance Explained (%)');
            % legend('Individual','Cumulative','Location','best');
            % xlim([0.5 numel(explained)+0.5]);
            % grid on;
            % hold off;

            trainResp_red = score(:,kStart:kComp);
            PCtunings = zeros(length(freqVal),kComp);
            PCvariation = zeros(length(freqVal),kComp);
            for nf = 1 : length(freqVal)
                PCtunings(nf,kStart:kComp) = mean(trainResp_red(trainFreq == freqVal(nf),:),1);
                PCvariation(nf,kStart:kComp) = std(trainResp_red(trainFreq == freqVal(nf),:),1);
            end
            
            mdl = fitcdiscr(trainResp_red,trainFreq,'DiscrimType','linear','Gamma',1);
            predFreq = predict(mdl,testResp*coeff(:,kStart:kComp));
    
            % collect output
            PCTall(nj,nt,nd,nc) = mean(predFreq == testFreq); % accuracy
            MDLall{nj,nt,nd,nc} = mdl; % model parameters
            PCFall{nj,nt,nd,nc} = cat(3, PCtunings, PCvariation); % accuracy
        end 
    end
    fprintf(newline)
end

end



%% PLOT CUMULATIVE VARIANCE EXPLAINED
figure('Color',[1 1 1],'Theme','Light'); hold on;
CVEall = squeeze(CVEall);
for nc = 1 : 3
    locCVmean = mean(CVEall(:,:,nc),1);
    locCVmean(locCVmean == 0) = [];
    locCVci = zeros(length(locCVmean),2);
    for ni = 1 : length(locCVmean)
        locCVci(ni,:) = mkCI(squeeze(CVEall(:,ni,nc)))';
    end
    % xVals = 1:length(locCVmean);
    % xLooped = [xVals flip(xVals)]; 
    % yLooped = cat(1,locCVci(:,1),flip(locCVci(:,2)));
    % fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.3,'EdgeColor','none');
    plot(nCellDraws(1:length(locCVmean)),locCVmean,'LineStyle','-','Marker','.','LineWidth',1,'Color',colors.(cellType(nc)));
end
set(gca,'XLim',[13 500],'YLim',[30 100],'YTick',10:10:100,'XTick',[15:60:435 500],'FontName','Arial');
ylabel('Variance explained (PC 1-15)','FontName','Arial');
xlabel('Number of neurons used','FontName','Arial');
title('Percentage of explained variance from 15 PCs','FontName','Arial');


%% PLOT PC SELECTIVITY STRENGTH BY NUMBER OF NEURONS
locSels = zeros(length(cellType),size(PCFall,3),size(PCFall,1)); % 3 x 33 x 200
np = 4; % target PC number

figure('Color',[1 1 1],'Theme','Light');
tiledlayout(2,2);
for np = 1 : 4
    nexttile(np); hold on;
    for nc = 1 : length(cellType) % LOOP CELL TYPE
        locCIs = zeros(2,size(PCFall,3));
        for nb = 1 : size(PCFall,3) % LOOP NEURON NUMBER BIN
            for ni = 1 : size(PCFall,1) % LOOP MODEL ITERATION
                PCtunings = PCFall{ni,1,nb,nc};
                if isempty(PCtunings)
                    break
                end
                locSels(nc,nb,ni) = quickSelectivity(PCtunings(:,np,1));
            end
            locCIs(:,nb) = mkCI(squeeze(locSels(nc,nb,:)));
        end
        locMeanVals = mean(squeeze(locSels(nc,:,:)),2);
        locIndx = locMeanVals ~= 0;
        xVals = nCellDraws(locIndx);
        xLooped = cat(2,xVals,flip(xVals));
        yLooped = cat(2,locCIs(1,locIndx),flip(locCIs(2,locIndx)));
        fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.3,'EdgeColor','none');
        plot(xVals,locMeanVals(locIndx),'LineStyle','-','Marker','.','LineWidth',1,'Color',colors.(cellType(nc)));
    end
    %set(gca,'XLim',[13 500],'YLim',[0.2 0.5],'XTick',[15:60:435 500],'FontName','Arial');
    set(gca,'XLim',[13 500],'XTick',[15:60:435 500],'FontName','Arial');
    ylabel('Frequency selectivity strength','FontName','Arial');
    xlabel('Number of neurons used','FontName','Arial');
    title(strcat("Average selectivity strength of PC",num2str(np)),'FontName','Arial');
    hold off;
end



%% PLOT PC VARIATION PER RESPONSE
locSTDs = zeros(length(cellType),size(PCFall,3),size(PCFall,1)); % 3 x 33 x 200
np = 4; % target PC number

figure('Color',[1 1 1],'Theme','Light');
tiledlayout(1,4);
for np = 1 : 4
    nexttile(np); hold on;
    for nc = 1 : length(cellType) % LOOP CELL TYPE
        locCIs = zeros(2,size(PCFall,3));
        for nb = 1 : size(PCFall,3) % LOOP NEURON NUMBER BIN
            for ni = 1 : size(PCFall,1) % LOOP MODEL ITERATION
                PCtunings = PCFall{ni,1,nb,nc};
                if isempty(PCtunings)
                    break
                end
                locSTDs(nc,nb,ni) = mean(PCtunings(:,np,2));
            end
            locCIs(:,nb) = mkCI(squeeze(locSTDs(nc,nb,:)));
        end
        locMeanVals = mean(squeeze(locSTDs(nc,:,:)),2);
        locIndx = locMeanVals ~= 0;
        xVals = nCellDraws(locIndx);
        xLooped = cat(2,xVals,flip(xVals));
        yLooped = cat(2,locCIs(1,locIndx),flip(locCIs(2,locIndx)));
        fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.3,'EdgeColor','none');
        plot(xVals,locMeanVals(locIndx),'LineStyle','-','Marker','.','LineWidth',1,'Color',colors.(cellType(nc)));
    end
    
    set(gca,'XLim',[13 500],'XTick',[15:60:435 500],'FontName','Arial');
    ylabel('Frequency selectivity strength','FontName','Arial');
    xlabel('Number of neurons used','FontName','Arial');
    title(strcat("Average selectivity strength of PC",num2str(np)),'FontName','Arial');
    hold off;
end

%% PLOT: PERFORMANCE BY CELL TYPE

% KEY: 
% [nj] Model Iterations     ... [scalar, 200]
% [nt] Selectivity Types    ... ["selective","random"]
% [nd] Number of Cell Draws ... [various lengths, 5-8 total]
% [nc] Cell Types           ... ["PV","SOM","EXC"]

figure('Color',[1 1 1],'Theme','Light'); hold on;
selType = ["selective","any"]; 
maxDraw = 500;
neuStepSize = 15; % num neurons to add to model each step
is.any = true(length(is.selective),1);

% fake data for correct legend order
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');
plot([-1 -2],[-1 -1],'k-');  % Freq Sel
plot([-1 -2],[-1 -1],'k--'); % Random Draws
plot([-1 -2],[-1 -1],'k:'); % Chance

% PLOT DATA FROM PCTall2 (FREQ. SEL NEURONS BINNED)
for nt = 1 : length(selType)
    for nc = 1 : length(cellType)
        if strcmp(cellType(nc),'EXC')
            maxNeuAvail = maxDraw;
        else % PV | SOM
            maxNeuAvail = sum(is.(cellType(nc)) & is.(selType(nt)));
            if maxNeuAvail > maxDraw
                maxNeuAvail = maxDraw;
            end
        end
        nCellDraws = 15:neuStepSize:maxNeuAvail;
        target = cellType(nc);
        nBins = 1:length(nCellDraws);
        yMean = mean(squeeze(PCTall(:,nt,nBins,nc)),1).*100;
        % ySems = sem(squeeze(PCTall2(targetAmp,:,nt,nBins,nc)),1).*100;
        yCIs = zeros(length(nBins),2); % confidence intervals (95%)
        for nd = 1 : length(nBins)
            yCIs(nd,:) = mkCI(squeeze(PCTall(:,nt,nd,nc))).*100;
        end
        xLooped = [nCellDraws flip(nCellDraws)]; 
        % yLooped = [yMean+ySems flip(yMean-ySems)]; 
        yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
        fill(xLooped,yLooped,colors.(target),'FaceAlpha',0.3,'EdgeColor','none');
        if nt == 1
            plot(nCellDraws,yMean,'-','LineWidth',1,'Color',colors.(target));
        elseif nt == 2
            plot(nCellDraws,yMean,'--','LineWidth',1,'Color',colors.(target));
        end

    end
end

plot([15 maxDraw],[20 20],'k:'); % chance performance level
ylabel('Frequency decoding accuracy (%)','FontName','Arial');
xlabel('Number of neurons used','FontName','Arial');
title('Neuron source improves model performance','FontName','Arial');
subtitle('All models use predictors: PC 1-15','FontName','Arial');
set(gca,'XLim',[13 500],'YLim',[18 100],'YTick',10:10:100,'XTick',[15:60:435 500]);
% legend({"Freq. Mod","All neurons","EXC","PV","SOM"},'Location','Northwest');
legend({"PV","SOM","EXC","Freq.selective","All available","Chance"},'Location','SouthEast',...
    'box','off','NumColumns',2);



%% PLOT: PERFORMANCE BY CELL TYPE
% PC2-15 (PV & EXC), PC3-15 (SOM)

% -----------------------------------------------
% cd('C:\Users\skich\Desktop\WORK');
% load("modelAccuracy_dataset1_PC2-15_PC3-15.mat");
% -----------------------------------------------

cellType = ["PV","SOM","EXC"];
selType = "selective";
figure('Color',[1 1 1],'Theme','Light'); hold on;

% fake data for correct legend order
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');
plot([-1 -2],[-1 -1],'k-');  % Inter & Dual
plot([-1 -2],[-1 -1],'k:'); % Freq & Amp

% PLOT DATA FROM PCTall2 (FREQ. SEL NEURONS BINNED)
for nt = 1 : length(selType)
    for nc = 1 : length(cellType)
        if strcmp(cellType(nc),'EXC')
            maxNeuAvail = maxDraw;
        else % PV | SOM
            maxNeuAvail = sum(is.(cellType(nc)) & is.(selType(nt)));
            if maxNeuAvail > maxDraw
                maxNeuAvail = maxDraw;
            end
        end
        nCellDraws = 15:neuStepSize:maxNeuAvail;
        target = cellType(nc);
        nBins = 1:length(nCellDraws);
        yMean = mean(squeeze(PCTall(:,nBins,nc)),1).*100;
        % ySems = sem(squeeze(PCTall2(targetAmp,:,nt,nBins,nc)),1).*100;
        yCIs = zeros(length(nBins),2); % confidence intervals (95%)
        for nd = 1 : length(nBins)
            yCIs(nd,:) = mkCI(squeeze(PCTall(:,nd,nc))).*100;
        end
        xLooped = [nCellDraws flip(nCellDraws)]; 
        % yLooped = [yMean+ySems flip(yMean-ySems)]; 
        yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
        fill(xLooped,yLooped,colors.(target),'FaceAlpha',0.3,'EdgeColor','none');
        plot(nCellDraws,yMean,'-','LineWidth',1,'Color',colors.(target));
    end
end


plot([15 maxDraw],[20 20],'k:'); % chance performance level
ylabel('Frequency decoding accuracy (%)','FontName','Arial');
xlabel('Neurons used','FontName','Arial');
title("Model performance with amp-correlated PC's removed",'FontName','Arial');
subtitle("Predictors: PC2-15 (PV & EXC), PC3-15 (SOM)",'FontName','Arial');
set(gca,'XLim',[13 500],'YLim',[18 100],'YTick',10:10:100,'XTick',15:60:500);
% legend({"Freq. Mod","All neurons","EXC","PV","SOM"},'Location','Northwest');
legend({"PV","SOM","EXC","Freq. selective","chance"},'numColumns',2,'Location','SouthEast','box','off');


% INSET BAR PLOT: AVG DIFFERENCE BETWEEN FREQ.SEL & ALL.AVAIL 
pvDiff = [];
somDiff = [];
excDiff = [];
for nn = 1 : size(PCTall,3)
    for nc = 1 : 3
        if all(PCTall(:,1,nn,nc) == 0) || all(PCTall(:,2,nn,nc) == 0)
            continue
        else
            % avg freq.sel value / all.avail value
            locDiff = (mean(PCTall(:,1,nn,nc)) - mean(PCTall(:,2,nn,nc)))*100; 
            switch nc
                case 1 % PV
                    pvDiff = cat(1,pvDiff,locDiff);
                case 2 % SOM
                    somDiff = cat(1,somDiff,locDiff);
                case 3 % EXC
                    excDiff = cat(1,excDiff,locDiff);
            end
        end
    end
end

figure('Color',[1 1 1]); hold on;
% PV BAR
bar(1,mean(pvDiff),0.8,'FaceColor',colors.PV,'EdgeColor',[0 0 0]);
mean(abs(mkCI(pvDiff)-mean(pvDiff)));
errorbar(1,mean(pvDiff),mean(abs(mkCI(pvDiff)-mean(pvDiff))),'vertical','LineWidth',1.2,'Color',[0 0 0]);
% PV BAR
bar(2,mean(somDiff),0.8,'FaceColor',colors.SOM,'EdgeColor',[0 0 0]);
mean(abs(mkCI(somDiff)-mean(somDiff)));
errorbar(2,mean(somDiff),mean(abs(mkCI(somDiff)-mean(somDiff))),'vertical','LineWidth',1.2,'Color',[0 0 0]);
% PV BAR
bar(3,mean(excDiff),0.8,'FaceColor',colors.EXC,'EdgeColor',[0 0 0]);
mean(abs(mkCI(excDiff)-mean(excDiff)));
errorbar(3,mean(excDiff),mean(abs(mkCI(excDiff)-mean(excDiff))),'vertical','LineWidth',1.2,'Color',[0 0 0]);
% LABELS
ylabel('\Delta Accuracy (%)','FontName','Arial');
title({'Frequency-selective','improvement'},'FontName','Arial');
set(gca,'XLim',[0.4 3.6],'XTick',[]);


%% PLOT: PERFORMANCE BY CELL TYPE
% PC2-15 (PV & EXC), PC3-15 (SOM)
% for gamma 0/1 data

cellType = ["PV","SOM","EXC"];
selType = "selective";
figure('Color',[1 1 1],'Theme','Light'); hold on;

% fake data for correct legend order
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');
plot([-1 -2],[-1 -1],'k-');  % Inter & Dual
plot([-1 -2],[-1 -1],'k:'); % Freq & Amp

% PLOT DATA FROM PCTall2 (FREQ. SEL NEURONS BINNED)
for nt = 1 : length(selType)
    for nc = 1 : length(cellType)
        if strcmp(cellType(nc),'EXC')
            maxNeuAvail = maxDraw;
        else % PV | SOM
            maxNeuAvail = sum(is.(cellType(nc)) & is.(selType(nt)));
            if maxNeuAvail > maxDraw
                maxNeuAvail = maxDraw;
            end
        end
        nCellDraws = 15:neuStepSize:maxNeuAvail;
        target = cellType(nc);
        nBins = 1:length(nCellDraws);
        yMean = mean(squeeze(PCTall(:,nt,nBins,nc)),1).*100;
        % ySems = sem(squeeze(PCTall2(targetAmp,:,nt,nBins,nc)),1).*100;
        yCIs = zeros(length(nBins),2); % confidence intervals (95%)
        for nd = 1 : length(nBins)
            yCIs(nd,:) = mkCI(squeeze(PCTall(:,nt,nd,nc))).*100;
        end
        xLooped = [nCellDraws flip(nCellDraws)]; 
        % yLooped = [yMean+ySems flip(yMean-ySems)]; 
        yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
        fill(xLooped,yLooped,colors.(target),'FaceAlpha',0.3,'EdgeColor','none');
        plot(nCellDraws,yMean,'-','LineWidth',1,'Color',colors.(target));
    end
end


plot([15 maxDraw],[20 20],'k:'); % chance performance level
ylabel('Frequency decoding accuracy (%)','FontName','Arial');
xlabel('Neurons used','FontName','Arial');
title("Model performance with amp-correlated PC's removed",'FontName','Arial');
subtitle("Predictors: PC2-15 (PV & EXC), PC3-15 (SOM)",'FontName','Arial');
set(gca,'XLim',[13 500],'YLim',[18 100],'YTick',10:10:100,'XTick',15:60:500);
% legend({"Freq. Mod","All neurons","EXC","PV","SOM"},'Location','Northwest');
legend({"PV","SOM","EXC","Freq. selective","chance"},'numColumns',2,'Location','SouthEast','box','off');




%% Other decoding approach
% pick max neuron draw (285), titrate gamma and delta

% perform population decoding vary 

nIter = 200; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
usePCA = true; % use PCA step prior to model training
kStart = 1; % first PC to include
kComp = 15; % total number of PCs to use in model training / testing
useAllNeu = false; % use all available neurons per cell type for modeling
cellType = ["PV","SOM","EXC"]; % order for plotting
% cellType = "SOM"; % order for plotting
is.any = true(length(is.selective),1);
selType = ["selective","any"]; 
% selType = ["selective"]; 
neuStepSize = 15; % num neurons to add to model each step
nCellDraws = 280;
gammaSteps = 0:0.2:1;
deltaSteps = 0:0.2:1;

% CREATE OUTPUT ARRAYS
PCTall = zeros(nIter,length(selType),length(cellType),length(gammaSteps),length(deltaSteps)); % (nj,nt,nc,ng,nd)
MDLall = cell(nIter,length(selType),length(cellType),length(gammaSteps),length(deltaSteps)); % (nj,nt,nc,ng,nd)
CVEall = zeros(nIter,length(selType),length(cellType),length(gammaSteps),length(deltaSteps));
PCFall = cell(nIter,length(selType),length(cellType),length(gammaSteps),length(deltaSteps)); % (nj,nt,nc,ng,nd)

% % Create output structure for storing PCA information
% PC = struct();
% for nc = 1 : length(cellType)
%     PC.(cellType(nc)) = struct();
% end


for nc = 1:length(cellType) %nt = 1 : length(selType)
disp(cellType(nc));
% disp(selType(nt));

if nc == 2
    kStart = 3;
else
    kStart = 2;
end

for nt = 1 : length(selType) %nd = 1 : length(nCellDraws)

    % nCellDraws = kComp : neuStepSize : maxNeuAvail;
    trainTrial = round(minTrialNum*trainPortion);
    testTrial = round(minTrialNum*(1-trainPortion));
    trainFreq = repmat(freqVal,trainTrial,1);
    testFreq = repmat(freqVal,testTrial,1);
    
    for nd = 1 : length(nCellDraws) %nc = 1:length(cellType)
        %fprintf(strcat(" ",cellType(nc)))
        fprintf(strcat(num2str(nCellDraws(nd))," "));
        nCellDraw = nCellDraws(nd);
        idx = find(Tc.identity == cellType(nc) & is.(selType(nt)));
        
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
            [coeff,score,latent] = pca(trainResp);
            % PC.(cellType(nc)).coeff{ni} = coeff;
            % PC.(cellType(nc)).score{ni} = score;
            % PC.(cellType(nc)).latent{ni} = latent;

            % Percent variance explained
            explained = 100 * latent / sum(latent);
            CVEall(nj,nt,nc,1,1) = sum(explained(1:15));
            
            % % Bar plot + cumulative line
            % figure;
            % bar(explained,'FaceColor',[0.6 0.6 0.9]); hold on;
            % plot(cumsum(explained),'r-o','LineWidth',1.5);
            % xlabel('Principal Component');
            % ylabel('Variance Explained (%)');
            % legend('Individual','Cumulative','Location','best');
            % xlim([0.5 numel(explained)+0.5]);
            % grid on;
            % hold off;

            trainResp_red = score(:,kStart:kComp);
            PCtunings = zeros(length(freqVal),kComp);
            PCvariation = zeros(length(freqVal),kComp);
            for nf = 1 : length(freqVal)
                PCtunings(nf,kStart:kComp) = mean(trainResp_red(trainFreq == freqVal(nf),:),1);
                PCvariation(nf,kStart:kComp) = std(trainResp_red(trainFreq == freqVal(nf),:),1);
            end
            
            for ng = 1 : length(gammaSteps)
                for ne = 1 : length(deltaSteps)
                    mdl = fitcdiscr(trainResp_red,trainFreq,'DiscrimType','linear','Gamma',gammaSteps(ng));
                    predFreq = predict(mdl,testResp*coeff(:,kStart:kComp));
                    % collect output
                    PCTall(nj,nt,nc,ng,ne) = mean(predFreq == testFreq); % accuracy
                    MDLall{nj,nt,nc,ng,ne} = mdl; % model parameters
                    PCFall{nj,nt,nc,ng,ne} = cat(3, PCtunings, PCvariation); % accuracy
                end
            end
        end 
    end
    fprintf(newline)
end
end

% per by gamma
figure('Color',[1 1 1]); hold on;
xSP = [-0.1,0,0.1];
for nc = 1 : length(cellType)
    allCIs =  zeros(2,length(gammaSteps));
    for ng = 1 : length(gammaSteps)
        allCIs(:,ng) = mkCI(squeeze(PCTall(:,1,nc,ng,1)));
        plot([ng ng]+xSP(nc), allCIs(:,ng)','-','Color',colors.(cellType(nc)));
        plot(ng+xSP(nc),mean(squeeze(PCTall(:,1,nc,ng,1))),'.','Color',colors.(cellType(nc)));
    end
end
set(gca,'XTick',1:6,'XTickLabels',0:0.2:1);
ylabel('Frequency decoding accuracy (%)','FontName','Arial');
xlabel('Gamma value','FontName','Arial');
title("Model performance by gamma",'FontName','Arial');


% per by delta
figure('Color',[1 1 1]); hold on;
xSP = [-0.1,0,0.1];
for nc = 1 : length(cellType)
    allCIs =  zeros(2,length(deltaSteps));
    for ng = 1 : length(deltaSteps)
        allCIs(:,ng) = mkCI(squeeze(PCTall(:,1,nc,1,ng)));
        plot([ng ng]+xSP(nc), allCIs(:,ng)','-','Color',colors.(cellType(nc)));
        plot(ng+xSP(nc),mean(squeeze(PCTall(:,1,nc,1,ng))),'.','Color',colors.(cellType(nc)));
    end
end
set(gca,'XTick',1:6,'XTickLabels',0:0.2:1);
ylabel('Frequency decoding accuracy (%)','FontName','Arial');
xlabel('Delta value','FontName','Arial');
title("Model performance by gamma",'FontName','Arial');