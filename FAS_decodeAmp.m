% FAS_decodeAmp
%
%
%

%% LOAD IN DATA

cd('C:\Users\skich\Desktop\WORK')
load("modelAccuracy_dataset2_ampDecode_selVsRandom.mat"); % FREQ.SEL vs ALL.NEURONS vs RAND.SHUFFLE

load("modelAccuracy_dataset2_ampDecode_allSelTypes.mat");

is.any = true(size(is.allMod));
cd('C:\Users\skich\Box\Tactile_Synchrony\Manuscript\Figures\Fig6')

%% perform population decoding

nIter = 200; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
targetFreq = 2; % 360, 525, 760, 1100, 1600 Hz
usePCA = true; % use PCA step prior to model training
kComp = 15; %number of PCs to use in model training / testing
useFreqEffNeu = true; % use neurons with significant anova freq effects
useAllNeu = false; % use all available neurons for modeling
includeErrBar = true; %include error bars in plots
cellType = ["PV","SOM","EXC"]; % order for plotting
% cellType = "EXC";
% selType = ["FreqMod","AmpMod","DualMod","InterMod"]; 

% % AMP SEL VS RANDOM DRAWS
is.allMod = (is.AmpMod | is.DualMod | is.InterMod) & ~is.trained;
% vis.any = true(size(is.allMod));
% selType = ["allMod","any"];
selType = "allMod";

neuStepSize = 5; % num neurons to add to model each step
numCellDraws = 1; % over-estimation to poulate arrays 60 (amp.sel vs. random) or 22 (4 sel types)

if length(selType) <= 2
    maxDraw = 160; % 2 SelTypes: FreqVsRandom (300), 
elseif length(selType) == 4 
    maxDraw = 120; % 4 SelTypes: Freq,Amp,Dual,Inter (120)
end

% FIGURE CREATION
% figure('Color',[1 1 1],'Theme','Light'); hold on;
%tiledlayout(2,2);

% CREATE OUTPUT ARRAYS
% Performance of Cell Types (all)
% PCTall = []; %for populating with accuracy measurements
% PCTall = zeros(nIter,length(selType),length(nCellDraws),length(cellType)); % (nj,nt,nd,nc)
PCTall = zeros(nFreq,nIter,length(selType),numCellDraws,length(cellType)); % (na,nj,nt,nd,nc)
% Models (all)
% MDLall = cell(nIter,length(selType),length(nCellDraws),length(cellType)); % (nj,nt,nd,nc)
MDLall = cell(nFreq,nIter,length(selType),numCellDraws,length(cellType)); % (na,nj,nt,nd,nc)

predAmpAll = cell(nFreq,nIter,length(selType),numCellDraws,length(cellType));
testAmpAll = cell(nFreq,nIter,length(selType),numCellDraws,length(cellType));


for nf = 1 : nFreq %num freqs
disp(strcat("FREQ VAL: ",num2str(freqVal(nf))));

for nc = 1:length(cellType) %nt = 1 : length(selType)
disp(cellType(nc));
% disp(selType(nt));
% determine number of neurons available
% PVnum = sum(is.PV & is.(selType(nt)));
% SOMnum = sum(is.SOM & is.(selType(nt)));
% maxNeuAvail = min([PVnum SOMnum]);
% nCellDraws = 15:5:maxNeuAvail;

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
    if numCellDraws == 1
        nCellDraws = maxNeuAvail;
    else
        nCellDraws = kComp:neuStepSize:maxNeuAvail;
    end
    trainTrial = round(minTrialNum*trainPortion);
    testTrial = round(minTrialNum*(1-trainPortion));
    trainFreq = repmat(freqVal,trainTrial,1);
    testFreq = repmat(freqVal,testTrial,1);
    trainAmp = repmat(ampValAll,trainTrial,1);
    testAmp = repmat(ampValAll,testTrial,1);
    
    tic;
    for nd = 1 : length(nCellDraws) %nc = 1:length(cellType)
        %fprintf(strcat(" ",cellType(nc)))
        fprintf(strcat(num2str(nCellDraws(nd))," "));
        nCellDraw = nCellDraws(nd);

        % select cell type
        if useFreqEffNeu == true
            if multiAmp == true
                idx = find(Tc.identity == cellType(nc) & is.(selType(nt)));
            elseif multiAmp == false
                idx = find(Tc.identity == cellType(nc) & Tc.selectivityPVal < 0.05);
            end
        elseif useFreqEffNeu == false
            idx = find(Tc.identity == cellType(nc));
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
                          
                for nl = 1:nAmp
                    % randomly draw trials at each amplitude, no replacement
                    if multiAmp
                        trialNum = find(cellFreq == freqVal(nf) & cellAmp == ampValAll(nl));
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
            trainResp = reshape(trainResp, nAmp*trainTrial,nCellDraw);
            testResp = reshape(testResp, nAmp*testTrial,nCellDraw);
    
            % Z-Scoring 
            for nz = 1 : nCellDraw
                mu = mean(trainResp(:,nz));
                sigma = std(trainResp(:,nz));
                trainResp(:,nz) = (trainResp(:,nz) - mu) ./ sigma;
                testResp(:,nz) = (testResp(:,nz) - mu) ./ sigma;
            end
    
            % Training
            % WITH PCA 
            [coeff,score] = pca(trainResp);
            % kComp = floor(trainTrial*nFreq*0.5);
            trainResp_red = score(:,1:kComp);
            mdl = fitcdiscr(trainResp_red,trainAmp,'DiscrimType','linear');
            % mdl = fitcdiscr(trainResp_red,trainAmp(randperm(length(trainAmp))),'DiscrimType','linear'); %train model (RAND MIXED DATA)
            predAmp = predict(mdl,testResp*coeff(:,1:kComp));
    
            % collect output
            PCTall(nf,nj,nt,nd,nc) = mean(predAmp == testAmp); %accuracy
            predAmpAll{nf,nj,nt,nd,nc} = predAmp;
            testAmpAll{nf,nj,nt,nd,nc} = testAmp;
            MDLall{nf,nj,nt,nd,nc} = mdl;
        end 
    end
    fprintf(newline)
end

end

end % nf (numFreqs)



%% PLOT: PERFORMANCE BY CELL TYPE AND SEL TYPE (AMP.MOD vs. RANDOM)
% ONE FREQ, EACH SELECTIVITY TYPE SEPARATE

targetFreq = 5;
maxDraw = 300;
selType = ["allMod","any","allMod"];
figure('Color',[1 1 1],'Theme','Light'); hold on;
% plot fake data for correct legend
fill([1 2 2 1],[1.1 1.1 .9 .9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
fill([1 2 2 1],[1.1 1.1 .9 .9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
fill([1 2 2 1],[1.1 1.1 .9 .9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');
plot([1 2],[1 1],'k-');
plot([1 2],[1 1],'k--');
%plot([1 2],[1 1],'k:');
neuStepSize = 5;
nCellDraws = 15:neuStepSize:maxDraw;

for nt = 1 : size(PCTall,3)
    for nc = 1 : length(cellType)
        useBins = ~all(squeeze(PCTall(targetFreq,:,nt,:,nc))==0,1);
        target = cellType(nc);
        yMean = mean(squeeze(PCTall(targetFreq,:,nt,useBins,nc)),1).*100;
        ySems = sem(squeeze(PCTall(targetFreq,:,nt,useBins,nc)),1).*100;
        xLooped = [nCellDraws(useBins) flip(nCellDraws(useBins))]; 
        yLooped = [yMean+ySems flip(yMean-ySems)]; 
        fill(xLooped,yLooped,colors.(target),'FaceAlpha',0.3,'EdgeColor','none');
        switch nt
            case 1
                plot(nCellDraws(useBins),yMean,'-','LineWidth',1,'Color',colors.(target));
            case 2 
                plot(nCellDraws(useBins),yMean,'--','LineWidth',1,'Color',colors.(target));
            case 3
                plot(nCellDraws(useBins),yMean,':','LineWidth',1,'Color',colors.(target));
        end
    end
end

ylabel('model accuracy (%)','FontName','Arial');
xlabel('Neurons used','FontName','Arial');
title(strcat("Amplitude decoding performance at ",num2str(freqVal(targetFreq))," hz"),'FontName','Arial');
% subtitle(strcat("Only trials at ",num2str(freqVal(targetFreq))," hz"),'FontName','Arial');
set(gca,'XLim',[13 maxDraw+2],'YLim',[20 75],'YTick',10:10:100,'XTick',15:30:maxDraw);
legend({"PV","SOM","EXC","Amp. selective","All available"},'NumColumns',2,...
    'box','off','Location','Northwest');



%% PERFORMANCE BY FREQUENCY (ALL-FREQ SEL vs. RANDOM DRAWS)

cd('D:\Tactile_Synchrony\2P-Data\Cleaned\datasets')
load("modelAccuracy_dataset2_ampDecode_selVsRandom.mat")

targetDrawBin = 28; % bin 28 is 150 neurons (15:5:300)
xVals = freqVal'; 
figure('Color',[1 1 1],'Theme','Light'); hold on;

% fake data for correct legend order
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');
% plot([-1 -2],[-1 -1],'k-');  % Inter & Dual
% plot([-1 -2],[-1 -1],'k--'); % Freq & Amp

% PLOT DATA FROM PCTall2 (FREQ. SEL NEURONS BINNED)
for nt = 1 : 1 % size(PCTall,3) % 2: sel types (freq-sel vs. random)
    for nc = 1 : size(PCTall,5) % 3: cell types (PV, SOM, EXC)
        target = cellType(nc);
        yMean = mean(squeeze(PCTall(:,:,nt,targetDrawBin,nc)),2).*100;
        yCIs = zeros(size(PCTall,1),2); % confidence intervals (95%)
        for na = 1 : size(PCTall,1) % 5: num Amps
            yCIs(na,:) = mkCI(squeeze(PCTall(na,:,nt,targetDrawBin,nc))).*100;
        end
        xLooped = [xVals flip(xVals)]; 
        % yLooped = [yMean+ySems flip(yMean-ySems)]; 
        yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
        fill(xLooped,yLooped,colors.(target),'FaceAlpha',0.3,'EdgeColor','none');
        switch nt
            case 1
                plot(xVals,yMean,'-','LineWidth',1,'Color',colors.(target));
            case 2 
                plot(xVals,yMean,'--','LineWidth',1,'Color',colors.(target));
        end
    end
end

ylabel('model accuracy (%)','FontName','Arial');
xlabel('Frequency (Hz)','FontName','Arial');
title('Amplitude decoding performance across frequencies','FontName','Arial');
subtitle('Fixed model size: 150 amplitude selective neurons','FontName','Arial');
set(gca,'XLim',[min(freqVal)-5 max(freqVal)+20],'YLim',[30 70],'YTick',10:5:100,...
    'XTick',xVals,'XTickLabelRotation',45,'XScale','log');
% legend({"PV","SOM","EXC","Amp. selective","All available"},'Location','Northeast',...
%     'box','off','NumColumns',2);
legend({"PV","SOM","EXC"},'Location','Northeast',...
    'box','off','NumColumns',1);




%% PLOT: PERFORMANCE BY CELL TYPE AND SEL TYPE
% TOP AMP, EACH SELECTIVITY TYPE SEPARATE

% KEY: 
% [nj] Model Iterations     ... [scalar, 200]
% [nt] Selectivity Types    ... ["FreqMod","AmpMod","DualMod","InterMod"]
% [nd] Number of Cell Draws ... [various lengths, 5-8 total]
% [nc] Cell Types           ... ["EXC","PV","SOM"]

targetFreq = 2;
maxDraw = 120; 
neuStepSize = 5; 
selType = ["FreqMod","AmpMod","DualMod","InterMod"]; 
figure('Color',[1 1 1],'Theme','Light'); hold on;
% plot fake data for correct legend
plot([1 2],[1 1],'k-'); plot([1 2],[1 1],'k-.'); % Inter & Dual
plot([1 2],[1 1],'k--'); plot([1 2],[1 1],'k:'); % Freq & Amp
plot([15 120], [20 20],'r--'); %chance
fill([1 2 2 1],[1.1 1.1 .9 .9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
fill([1 2 2 1],[1.1 1.1 .9 .9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
fill([1 2 2 1],[1.1 1.1 .9 .9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');

for nt = 1 : length(selType)
    for nc = 1 : length(cellType)
        if strcmp(cellType(nc),'EXC')
            maxNeuAvail = maxDraw;
        else % PV | SOM
            maxNeuAvail = sum(is.(cellType(nc)) & is.(selType(nt)));
        end
        nCellDraws = 15:neuStepSize:maxNeuAvail;
        target = cellType(nc);
        nBins = 1:length(nCellDraws);
        yMean = mean(squeeze(PCTall(targetFreq,:,nt,nBins,nc)),1).*100;
        ySems = sem(squeeze(PCTall(targetFreq,:,nt,nBins,nc)),1).*100;
        xLooped = [nCellDraws flip(nCellDraws)]; 
        yLooped = [yMean+ySems flip(yMean-ySems)]; 
        fill(xLooped,yLooped,colors.(target),'FaceAlpha',0.3,'EdgeColor','none');
        switch selType(nt)
            case "FreqMod"
                plot(nCellDraws,yMean,'--','LineWidth',1,'Color',colors.(target));
            case "AmpMod"
                plot(nCellDraws,yMean,':','LineWidth',1,'Color',colors.(target));
            case "DualMod"
                plot(nCellDraws,yMean,'-.','LineWidth',1,'Color',colors.(target));
            case "InterMod" 
                plot(nCellDraws,yMean,'-','LineWidth',1,'Color',colors.(target));
            case "allMod" 
                plot(nCellDraws,yMean,'-','LineWidth',1,'Color',colors.(target));    
        end
    end
end

ylabel('model accuracy (%)','FontName','Arial');
xlabel('Number of neurons used','FontName','Arial');
title('Amplitude decoding accuracy by effect type (525 Hz)','FontName','Arial');
% subtitle(strcat("Only trials at ",num2str(freqVal(targetFreq))," hz"),'FontName','Arial');
set(gca,'XLim',[15 maxDraw],'YLim',[15 70],'YTick',20:10:70,'XTick',15:15:120);

% % ADDITIONAL FIG FOR LEGEND
% figure('theme','light'); hold on;
% plot([1 2],[1 1],'k-'); plot([1 2],[1 1],'k-.'); % Inter & Dual
% plot([1 2],[1 1],'k--'); plot([1 2],[1 1],'k:'); % Freq & Amp
% fill([1 2 2 1],[1.1 1.1 .9 .9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
% fill([1 2 2 1],[1.1 1.1 .9 .9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
% fill([1 2 2 1],[1.1 1.1 .9 .9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');
% plot([1 2],[1 1],'r--'); % chance performance level
% set(gca,'XLim',[1 5],'YLim',[1 5]);
% legend({"F*A","F+A","F","A","PV","SOM","EXC","chance"},'NumColumns',2,...
%    'Location','northwest','box','off');



%% PLOT: PERFORMANCE OF EXC AND SEL TYPE OVER 1000 NEURONS
% TOP AMP, EACH SELECTIVITY TYPE SEPARATE

% KEY: 
% [nj] Model Iterations     ... [scalar, 200]
% [nt] Selectivity Types    ... ["FreqMod","AmpMod","DualMod","InterMod"]
% [nd] Number of Cell Draws ... [various lengths, 5-8 total]
% [nc] Cell Types           ... ["EXC","PV","SOM"]

targetFreq = 2;
maxDraw = 1000; % 120 or 1000
neuStepSize = 15; % 5 or 15
selType = ["FreqMod","AmpMod","DualMod","InterMod"]; 
figure('Color',[1 1 1],'Theme','Light'); hold on;
% plot fake data for correct legend
plot([1 2],[1 1],'k-'); plot([1 2],[1 1],'k-.'); % Inter & Dual
plot([1 2],[1 1],'k--'); plot([1 2],[1 1],'k:'); % Freq & Amp
plot([15 maxDraw], [20 20],'r--'); %chance
%fill([1 2 2 1],[1.1 1.1 .9 .9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
%fill([1 2 2 1],[1.1 1.1 .9 .9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
%fill([1 2 2 1],[1.1 1.1 .9 .9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');

for nt = 1 : length(selType)
    for nc = 1 : length(cellType)
        if strcmp(cellType(nc),'EXC')
            maxNeuAvail = maxDraw;
        else % PV | SOM
            maxNeuAvail = sum(is.(cellType(nc)) & is.(selType(nt)));
        end
        nCellDraws = 15:neuStepSize:maxNeuAvail;
        target = cellType(nc);
        nBins = 1:length(nCellDraws);
        yMean = mean(squeeze(PCTall(targetFreq,:,nt,nBins,nc)),1).*100;
        ySems = sem(squeeze(PCTall(targetFreq,:,nt,nBins,nc)),1).*100;
        xLooped = [nCellDraws flip(nCellDraws)]; 
        yLooped = [yMean+ySems flip(yMean-ySems)]; 
        fill(xLooped,yLooped,colors.(target),'FaceAlpha',0.3,'EdgeColor','none');
        switch selType(nt)
            case "FreqMod"
                plot(nCellDraws,yMean,'--','LineWidth',1,'Color',colors.(target));
            case "AmpMod"
                plot(nCellDraws,yMean,':','LineWidth',1,'Color',colors.(target));
            case "DualMod"
                plot(nCellDraws,yMean,'-.','LineWidth',1,'Color',colors.(target));
            case "InterMod" 
                plot(nCellDraws,yMean,'-','LineWidth',1,'Color',colors.(target));
            case "allMod" 
                plot(nCellDraws,yMean,'-','LineWidth',1,'Color',colors.(target));    
        end
    end
end

ylabel('model accuracy (%)','FontName','Arial');
xlabel('Number of neurons used','FontName','Arial');
title('Amplitude decoding accuracy by effect type (525 Hz)','FontName','Arial');
subtitle("EXC neurons only",'FontName','Arial');
set(gca,'XLim',[15 maxDraw],'YLim',[15 90],'YTick',20:10:90,'XTick',15:75:1000);

% % ADDITIONAL FIG FOR LEGEND
% figure('theme','light'); hold on;
% plot([1 2],[1 1],'k-'); plot([1 2],[1 1],'k-.'); % Inter & Dual
% plot([1 2],[1 1],'k--'); plot([1 2],[1 1],'k:'); % Freq & Amp
% fill([1 2 2 1],[1.1 1.1 .9 .9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
% fill([1 2 2 1],[1.1 1.1 .9 .9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
% fill([1 2 2 1],[1.1 1.1 .9 .9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');
% plot([1 2],[1 1],'r--'); % chance performance level
% set(gca,'XLim',[1 5],'YLim',[1 5]);
legend({"F*A","F+A","F","A","chance"},'NumColumns',2,...
    'Location','northwest','box','off');


%% ACCURACY AT N=40, EACH OF 4 SEL TYPES
targetFreq = 1; 
targetBin = 6; % 40 cell draws
xSpace = [1,2,3;6,7,8;11,12,13;16,17,18];
ntOrder = [2 1 3 4];
ntNames = {"Amplitude","Frequency","Dual Effect","Interaction"};
figure('theme','light','color',[1 1 1]); hold on;
plot(-1,-1,'.','MarkerSize',8,'Color',colors.PV);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.SOM);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.EXC);
for nt = 1 : 4
    for nc = 1 : length(cellType)
        yVals = PCTall(targetFreq,:,ntOrder(nt),targetBin,nc);
        yMean = mean(yVals)*100; % mean y value fonverted to percentage
        yCIs = (mkCI(yVals) - mean(yVals)).*100; % 95% confidence intervals centered at mean
        plot(xSpace(nt,nc),yMean,'.','MarkerSize',8,'Color',colors.(cellType(nc)));
        errorbar(xSpace(nt,nc),yMean,yCIs(1),yCIs(2),'Color',colors.(cellType(nc)));
    end
end

ylabel('model accuracy (%)','FontName','Arial');
xlabel('Selectivity Type','FontName','Arial');
title('Amplitude decoding accuracy by selectivity type','FontName','Arial');
subtitleText = strcat("40 neurons per condition, frequency: ",num2str(freqVal(targetFreq)),"Hz");
subtitle(subtitleText,'FontName','Arial');
set(gca,'YLim',[15 55],'XLim',[0 18],'XTick',xSpace(:,2),'XTickLabel',ntNames,'FontName','Arial');
legend({"PV","SOM","EXC"},'Location','Northwest','FontName','Arial','box','off');

%% ACCURACY AT N = 40, EACH OF 4 CELL TYPES, MAX & MIN AMP
targetFreq = [1 5]; % 360 & 1600 Hz
targetBin = 6; % 40 cell draws
xSpace = [1,2,3;7,8,9;13,14,15;19,20,21];
xGap = 0.5;
ntOrder = [1 2 3 4];
ntNames = {"Frequency","Amplitude","Dual Effect","Interaction"};
figure('theme','light','color',[1 1 1]); hold on;
plot(-1,-1,'.','MarkerSize',8,'Color',colors.PV);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.SOM);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.EXC);
plot(-1,-1,'o','MarkerSize',5,'Color',colors.EXC);
plot(-1,-1,'.','MarkerSize',11,'Color',colors.EXC);

for nc = 1 : length(cellType)
    yMeans = zeros(4,2);
    for nt = 1 : length(selType)
        for nf = 1 : length(targetFreq)
            yVals = PCTall(targetFreq(nf),:,ntOrder(nt),targetBin,nc);
            yMeans(nt,nf) = mean(yVals)*100; % mean y value fonverted to percentage
            yCIs = (mkCI(yVals) - mean(yVals)).*100; % 95% confidence intervals centered at mean
            if nf == 1
                plot(xSpace(nt,nc),yMeans(nt,nf),'o','MarkerSize',5,'Color',colors.(cellType(nc)));
                errorbar(xSpace(nt,nc),yMeans(nt,nf),yCIs(1),yCIs(2),'Color',colors.(cellType(nc)));
            elseif nf == 2
                plot(xSpace(nt,nc)+xGap,yMeans(nt,nf),'.','MarkerSize',11,'Color',colors.(cellType(nc)));
                errorbar(xSpace(nt,nc)+xGap,yMeans(nt,nf),yCIs(1),yCIs(2),'Color',colors.(cellType(nc)));
            end
        end
    end
end

ylabel('model accuracy (%)','FontName','Arial');
xlabel('Selectivity Type','FontName','Arial');
title('Amplitude decoding accuracy by selectivity type','FontName','Arial');
subtitle('40 neurons per condition','FontName','Arial');
set(gca,'YLim',[17 55],'XLim',[0 22],'XTick',xSpace(:,2)+(xGap/2),'XTickLabel',ntNames,'FontName','Arial');
legend({"PV","SOM","EXC","360 Hz","1600 Hz"},'Location','Northwest','FontName','Arial',...
    'box','off','NumColumns',2);


%% ACCURACY AT N=40, EACH OF 4 CELL TYPES, AVERAGE ACROSS FREQS
targetBin = 6; % 40 cell draw size
nFreqs = size(PCTall,1);
xSpace = [1,2,3;7,8,9;13,14,15;19,20,21];
xGap = 0.5;
ntOrder = [1 2 3 4];
ntNames = {"Frequency","Amplitude","Dual Effect","Interaction"};
figure('theme','light','color',[1 1 1]); hold on;
plot(-1,-1,'.','MarkerSize',8,'Color',colors.PV);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.SOM);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.EXC);


for nc = 1 : length(cellType)
    yMeans = zeros(4,nFreqs);
    for nt = 1 : length(ntNames)
        %for na = 1 : length(targetAmp)
            yVals=PCTall(:,:,ntOrder(nt),targetBin,nc);
            yMeans(nt,:) = mean(yVals,2).*100; % mean y value converted to percentage
            yMean = mean(yVals,'all')*100;
            yCI = mkCI(yVals(:)).*100;
            yCIs = zeros(size(yVals,1),2);
            for nf = 1 : nFreqs
                yCIs(nf,:) = (mkCI(yVals(nf,:)) - mean(yVals(nf,:))).*100; % 95% confidence intervals centered at mean
            end
            plot(xSpace(nt,nc),yMean,'.','MarkerSize',8,'Color',colors.(cellType(nc)));
            errorbar(xSpace(nt,nc),yMean,yCI(1)-yMean,yCI(2)-yMean,'Color',colors.(cellType(nc)));
        %end
    end
end


ylabel('model accuracy (%)','FontName','Arial');
xlabel('Selectivity Type','FontName','Arial');
title('Amplitude decoding performance of each selectivity type','FontName','Arial');
subtitle('40 neurons per condition, averaged across frequency','FontName','Arial');
set(gca,'YLim',[17 53],'XLim',[0 22],'XTick',xSpace(:,2)+(xGap/2),'XTickLabel',ntNames,'FontName','Arial');
legend({"PV","SOM","EXC"},'Location','Northwest','FontName','Arial',...
    'box','off');


%% ACCURACY AT N=40, EACH OF 4 SEL TYPES, AVG BETWEEN CELL TYPES
targetFreq = 2; 
targetBin = 6; % 40 cell draws
xSpace = [1,2,3,4,5; 8,9,10,11,12; 15,16,17,18,19; 22,23,24,25,26];
ntOrder = [2 1 3 4];
ntNames = {"Amplitude","Frequency","Dual Effect","Interaction"};

%cmap = colormap("parula");
%faceColors = cmap(round(linspace(1,256,5)),:); % matlab default
faceColors = mkcolors([0.8,0.8,0.8],[0 0 0],5);

figure('theme','light','color',[1 1 1]); hold on;
plot(-1,-1,'.','MarkerSize',8,'Color',faceColors(1,:));
plot(-1,-1,'.','MarkerSize',8,'Color',faceColors(2,:));
plot(-1,-1,'.','MarkerSize',8,'Color',faceColors(3,:));
plot(-1,-1,'.','MarkerSize',8,'Color',faceColors(4,:));
plot(-1,-1,'.','MarkerSize',8,'Color',faceColors(5,:));

for nf = 1 : 5
for nt = 1 : 4
    yVals = squeeze(PCTall(nf,:,ntOrder(nt),targetBin,nc));
    yMean = mean(yVals(:))*100; % mean y value fonverted to percentage
    yCIs = (mkCI(yVals(:)) - mean(yVals(:))).*100; % 95% confidence intervals centered at mean
    plot(xSpace(nt,nf),yMean,'.','MarkerSize',8,'Color',faceColors(nf,:));
    errorbar(xSpace(nt,nf),yMean,yCIs(1),yCIs(2),'Color',faceColors(nf,:));
end
end

ylabel('model accuracy (%)','FontName','Arial');
xlabel('Selectivity Type','FontName','Arial');
title('Amplitude decoding accuracy by selectivity type','FontName','Arial');
subtitleText = strcat("40 neurons per condition, frequency: ",num2str(freqVal(targetFreq)),"Hz");
subtitle(subtitleText,'FontName','Arial');
set(gca,'YLim',[15 55],'XLim',[0 27],'XTick',xSpace(:,2),'XTickLabel',ntNames,'FontName','Arial');
legend(string(freqVal),'Location','Northwest','FontName','Arial','box','off');


%% ACCURACY BY AMP 3x SEL TYPES (PSEUDO AMP XVAL & LINEAR FIT)
% KEY: 
% [nj] Model Iterations     ... [scalar, 200]
% [nt] Selectivity Types    ... ["FreqMod","AmpMod","DualMod","InterMod"]
% [nd] Number of Cell Draws ... [various lengths, 5-8 total]
% [nc] Cell Types           ... ["EXC","PV","SOM"]


targetBin = 6; % 8th bin = 50 neuron draws
xVals = freqVal'; % pseudo freqs
targetSels = [2 3 4]; % freq, dual, inter
figure('Color',[1 1 1],'Theme','Light'); hold on;
for nt = 1 : length(targetSels) % 3 sel types
    for nc = 1 : length(cellType)
        % MDL KEY:  nFreq (5), nIter (200), nSelTypes (3), nDrawBins (5-38), nCellTypes (3)
        yMean = mean(squeeze(PCTall(:,:,targetSels(nt),targetBin,nc))',1).*100;
        ySems = sem(squeeze(PCTall(:,:,targetSels(nt),targetBin,nc))',1).*100;
        xLooped = [xVals flip(xVals)]; 
        yLooped = [yMean+ySems flip(yMean-ySems)]; 
        fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.3,'EdgeColor','none');
        switch nt
            case 1 % freq
                plot(xVals,yMean,':','LineWidth',1,'Color',colors.(cellType(nc)));  
            case 2 % dual
                plot(xVals,yMean,'-.','LineWidth',1,'Color',colors.(cellType(nc)));  
            case 3 %inter
                plot(xVals,yMean,'-','LineWidth',1,'Color',colors.(cellType(nc)));  
        end
    end
end
% plot([xVals(1) xVals(end)],[20 20],'r--'); %chance performance line
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Frequency (Hz)','FontName','Arial');
titleTxt = "Amplitude decoding performance by frequency";
title(titleTxt,'FontName','Arial');
subtitle("From three selectivity types, 40 neurons per model")
set(gca,'XLim',[xVals(1)-10 xVals(end)+10],'YLim',[20 55],'YTick',20:10:60,...
    'XTick',xVals,'XTickLabel',freqVal);

%%  FITTING PLOT
f = struct(); gof = struct();
figure('Color',[1 1 1],'Theme','Light'); hold on;
% fake data for correct legend order
plot([-1 -2],[-1 -1],'k-'); plot([1 2],[1 1],'k-.'); % Inter & Dual
plot([-1 -2],[-1 -1],'k--'); % Freq & Amp
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');
set(gca,'XLim',[1 5],'YLim',[1 5]);
disp("-------------------------------------");

for nt = 1 : length(targetSels)
    for nc = 1 : length(cellType)
        yMeans.(cellType(nc)) = mean(squeeze(PCTall(:,:,targetSels(nt),targetBin,nc)),2);
        [locFit, locGof] = fit(xVals', yMeans.(cellType(nc)),'poly1');
        f.(cellType(nc)).(selType(targetSels(nt))) = locFit;
        gof.(cellType(nc)).(selType(targetSels(nt))) = locGof;
        switch nt
            case 1 % freq
                plot(xVals',locFit(xVals).*100,'--','Color',colors.(cellType(nc)));
            case 2 % dual
                plot(xVals',locFit(xVals).*100,'-.','Color',colors.(cellType(nc)));
            case 3 % inter
                plot(xVals',locFit(xVals).*100,'-','Color',colors.(cellType(nc)));
        end
        locStr = strcat(cellType(nc),"-",selType(targetSels(nt)),": rsq=",...
            num2str(gof.(cellType(nc)).(selType(targetSels(nt))).rsquare,3)," m=",...
            num2str(f.(cellType(nc)).(selType(targetSels(nt))).p1,3), " b=",...
            num2str(f.(cellType(nc)).(selType(targetSels(nt))).p2,3));
        disp(locStr);
    end
end
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Amplitude (\mum)','FontName','Arial');
titleTxt = "Frequency decoding performance by amplitude";
title(titleTxt,'FontName','Arial');
subtitle("Linear Fits");
set(gca,'XLim',[xVals(1) xVals(end)],'YLim',[25 70],'YTick',30:10:80,...
    'XTick',xVals,'XTickLabel',ampValAll);
legend({"Interaction","Both mod.","Freq. mod.","PV","SOM","EXC"},'NumColumns',2);
% legend({"","EXC","","PV","","SOM","fit","fit","fit"},'Location','NorthWest','NumColumns',2);

% % ADDITIONAL FIG FOR LEGEND
% figure('theme','light'); hold on;
% plot([-1 -2],[-1 -1],'k-'); plot([1 2],[1 1],'k-.'); % Inter & Dual
% plot([-1 -2],[-1 -1],'k--'); % Freq & Amp
% fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
% fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
% fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');
% set(gca,'XLim',[1 5],'YLim',[1 5]);
% legend({"Interaction","Both mod.","Freq. mod.","PV","SOM","EXC"},'NumColumns',2);


%% BIAS OVER/UNDER PREDICTION ANALYSIS

% cd('D:\Tactile_Synchrony\2P-Data\Cleaned\datasets');
% load('modelResults_dataset2_ampDecode_predictionAnalysis.mat'); % 'PCTall','MDLall','predFreqAll','testFreqAll'
% dimension IDs: nFreq, nIter, selType, 1, cellType

targSel = 1;
xSp = [-0.2 0 0.2]; % x-axis plot spacing
allDiffs = zeros(size(PCTall,1),size(PCTall,2),size(PCTall,3),size(PCTall,5)); % 5 x 200 x 4 x 3

figure('Color',[1 1 1]); hold on;
for nf = 1 : size(PCTall,1) % loop frequency (5)
    for nc = 1 : size(PCTall,5) % loop Cell Type (3)
        for ni = 1 : size(PCTall,2) % loop nIter (200)
            locPredVals = predAmpAll{nf,ni,targSel,1,nc}; % predicted amplitudes
            locTrueVals = testAmpAll{nf,ni,targSel,1,nc}; % true amplitudes
            locDiffs = locPredVals - locTrueVals; % differences (predicted - true)
            locDiffs(locDiffs == 0) = []; % remove correct values (0)
            allDiffs(nf,ni,targSel,nc) = mean(locDiffs);
        end
        locX = nf + xSp(nc); % x-position
        locY = mean(squeeze(allDiffs(nf,:,targSel,nc)),'omitnan'); % average difference
        locCI = mkCI(squeeze(allDiffs(nf,:,targSel,nc))); % 95 percent condifidence intervals
        errorbar(locX,locY,locY-locCI(1),'vertical','Color',colors.(cellType(nc)),'LineWidth',1);
        plot(locX,locY,'.','Color',colors.(cellType(nc)));
    end
end
set(gca,'XTick',1:5,'XTickLabels',freqValAll,'XLim',[0.6 5.4]);
xlabel('Frequency (Hz)','FontName','Arial');
ylabel({'Average prediction error (\mum)','(predicted - true)'},'FontName','Arial');
title('Amplitude decoding bias','FontName','Arial');

% EVALUATE WITH MODEL FITS AND PLOT
for nc = 1 : length(cellType)
    disp(cellType(nc));
    locYs = squeeze(allDiffs(:,:,targSel,nc)); % average difference values, amp by reps
    locY_vect = locYs(~isnan(locYs)); % remove NaN and vectorize
    locXs = repmat((1:5)',1,size(allDiffs,2)); % average difference values, amp by reps
    locX_vect = locXs(~isnan(locYs)); % remove NaN and vectorize

    % LINEAR FIT
    locMdl = fitlm(locX_vect,locY_vect); % linear fit
    AIC_linear = locMdl.ModelCriterion.AIC; % AIC
    fitSigLinear = locMdl.ModelFitVsNullModel.Pvalue < pValThresh; % significance of linear fit
    slope_Linear = locMdl.Coefficients.Estimate(2);
    rsq_Linear = locMdl.Rsquared.Adjusted;
    pval_Linear = locMdl.ModelFitVsNullModel.Pvalue;

    % EXP 1 FIT
    [fitobject,gof,output] = fit(locX_vect,locY_vect,'exp2');
    AIC_exp = output.numobs * log(gof.sse / output.numobs) + 2 * output.numparam; % AIC
    confBounds = confint(fitobject); % 95% confidence bounds of coefficients
    fitSigExp = all(confBounds(1,:).*confBounds(2,:) > 0); % check for same sign in coeff bounds

    % COMPARE MODELS
    xFine = linspace(min(locX_vect),max(locX_vect),100)';
    if fitSigLinear == true && fitSigExp == false
        % ONLY Linear model significant
        yFit = predict(locMdl,xFine);
        plot(xFine+xSp(nc),yFit,'--','Color',colors.(cellType(nc)));
    elseif fitSigLinear == false && fitSigExp == true 
        % ONLY Exponential model significant
        yFit = feval(fitobject, xFine); 
        yPI = predint(fitobject, xFine, 0.95); % 95% prediction intervals
        coefBounds = confint(fitobject, 0.95); % coeff. conf. bounds
        plot(xFine+xSp(nc),yFit,'--','Color',colors.(cellType(nc)));
    elseif fitSigLinear == true && fitSigExp == true 
        % BOTH models are significant, compare AIC
        if AIC_linear > AIC_exp 
            % Linear model has higher AIC
            yFit = predict(locMdl,xFine);
            plot(xFine+xSp(nc),yFit,'--','Color',colors.(cellType(nc)));
        elseif AIC_linear < AIC_exp 
            % Exponential model has higher AIC
            yFit = feval(fitobject, xFine); 
            plot(xFine+xSp(nc),yFit,'--','Color',colors.(cellType(nc)));
        else 
            % EQUAL AIC value (unlikely)
            error('models have equal AIC');
        end
    else 
        % NEITHER model is significant
        error('Both models are insignificant fits');
    end
    disp(strcat(cellType(nc),": (p-val. = ",num2str(pval_Linear,3),", r-sq. = ",...
        num2str(rsq_Linear,3),", slope = ",num2str(slope_Linear,3),")"));
end












