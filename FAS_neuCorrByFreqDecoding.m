
% FAS_neuCorrByFreqDecoding.m
%
% 

%% FOR DATASET 1

cd('D:\Tactile_Synchrony\2P-Data\Cleaned\datasets'); 

load("PCT_dataset1_workspace_20260616.mat"); % 'PCTall_baseline','PCTall_late','PCTallEXC_baseline','PCTallEXC_late'
load("PCT_dataset1_allAvailableNeurons.mat") % 'PCTall_baseline_allAvail','PCTallEXC_baseline_allAvail'

load("dataset1_neuCorrValsPerSess_normalRespWindows.mat"); % 'allCorrVals' & 'allPairNums'
allCorrValsBase = allCorrVals;
allPairNumsBase = allPairNums;

load("dataset1_neuCorrValsPerSess_lateRespWindows.mat"); % 'allCorrVals' & 'allPairNums'
allCorrValsLate = allCorrVals;
allPairNumsLate = allPairNums;

%% FOR DATASET 2
% cd('C:\Users\skich\Desktop\WORK');
% load('PCT_workspace_20260610.mat'); % PCTall_baseline, PCTall_late, PCTallExc_baseline, & PCTallExc_late
% 
% load('dataset2_neuCorrValsPerSess_lateRespWindows.mat'); % allCorrVals & allPairNums
% allCorrValsLate = allCorrVals;
% allPairNumsLate = allPairNums;
% 
% load('dataset2_neuCorrValsPerSess.mat'); % allCorrVals & allPairNums
% allCorrValsBase = allCorrVals;
% allPairNumsBase = allPairNums;
% 
% clear allCorrVals allPairNums 


%% 

% User Inputs
corrThreshDataset1 = 1; % within-sess threshold for extremely high corr values
corrThreshDataset2 = 0.3; % within-sess threshold for extremely high corr values

% FILTER OUT TRAINED SESSIONS
pairNames = {"EXC:EXC", "EXC:PV", "EXC:SOM", "PV:PV", "SOM:SOM"};
trainedSess = false(size(is.sessID,2),1);
for ns = 1 : size(is.sessID,2)
    locNeuList = find(is.sessID(:,ns));
    trainedSess(ns) = Tc.isTrained(locNeuList(1));
end

accMat = zeros(2,length(cellType)); % col.1 (acc mean, acc CI), col.2 (PV, SOM, EXC), Freq.Sel. Neurons
accMatAA = zeros(2,length(cellType)); % col.1 (acc mean, acc CI), col.2 (PV, SOM, EXC), All-Avail. Neurons
allAcc = cell(length(cellType),1);
allAccAA = cell(length(cellType),1);

% EXC:EXC BASELINE WINDOW: (0-1sec post-stim) 
PCTloc = squeeze(mean(PCTallEXC_baseline(:,~trainedSess,1),1))' .* 100; % frequency selective neurons
PCTloc(PCTloc == 0) = NaN;
PCTlocAA = squeeze(mean(PCTallEXC_baseline_allAvail(:,~trainedSess,1),1))' .* 100; %all-available neurons
PCTlocAA(PCTlocAA == 0) = NaN;

% SPLIT PAIRWISE CORRELATIONS UP BY SESSION & CELL TYPE (PV|SOM|EXC)
allPairTypes = cell(size(allCorrValsBase));
for ns = 1 : size(allPairNumsBase,1)
    allPairTypes{ns} = zeros(length(allPairNumsBase{ns}),1);
    for np = 1 : length(allPairNumsBase{ns})
        locPairID = Tc.identity(allPairNumsBase{ns}(np,:));
        if all(strcmp(locPairID,'EXC')) % EXC & EXC [1]
            allPairTypes{ns}(np) = 1;
        elseif any(strcmp(locPairID,'EXC')) && any(strcmp(locPairID,'PV')) % EXC & PV [2]
            allPairTypes{ns}(np) = 2;
        elseif any(strcmp(locPairID,'EXC')) && any(strcmp(locPairID,'SOM')) % EXC & SOM [3]
            allPairTypes{ns}(np) = 3;
        elseif all(strcmp(locPairID,'PV')) % PV & PV [4]
            allPairTypes{ns}(np) = 4;
        elseif all(strcmp(locPairID,'SOM')) % SOM & SOM [5]
            allPairTypes{ns}(np) = 5;
        end
    end
end

% PREPARE EXC:EXC PAIR CORRELATIONS / CELL NUMBERS
npt = 1; % EXC:EXC pairs only
locPairVals = NaN(size(allPairNumsBase,1),1);
locMeanCorr1 = NaN(size(allPairNumsBase,1),1);
locInhType = zeros(size(allPairNumsBase,1),1);
for ns = 1 : size(allPairNumsBase,1) % loop sessions
    if multiAmp == true
        locIdx = all(is.allMod(allPairNumsBase{ns}),2); % neuron pairs within-sess where both are selective
    elseif multiAmp == false
        locIdx = all(is.selective(allPairNumsBase{ns}),2); % neuron pairs within-sess where both are selective
    end
    locVals = allCorrValsBase{ns}(allPairTypes{ns} == npt & locIdx);
    locPairVals(ns) = mean(locVals);
    locMeanCorr1(ns) = mean(locVals); % corr values from selective pairs
    if any(allPairTypes{ns}==2) % PV present
        locInhType(ns) = 1;
    elseif any(allPairTypes{ns}==3) % SOM present
        locInhType(ns) = 2;
    end
end
locMeanCorr1 = locMeanCorr1(~trainedSess); % remove trained sessions

% APPLY THRESHOLDS TO REMOVE SESSION WITH HIGH NEURON CORRELATIONS 
% (per MGR suggestion)
if multiAmp == true
    abvThresh = locMeanCorr1 > corrThreshDataset2;
elseif multiAmp == false
    abvThresh = locMeanCorr1 > corrThreshDataset1;
end
locMeanCorr1(abvThresh) = [];
PCTloc1 = PCTloc(~abvThresh);
PCTloc1AA = PCTlocAA(~abvThresh);

% PLOTTING (EXC:EXC AT BASELINE)
figure('Color',[1 1 1]);
hold on
locInhType(abvThresh) = [];
locCorrPV = []; locAccPV = [];
locCorrSOM = []; locAccSOM = [];
for np = 1 : length(locMeanCorr1)
    if locInhType(np) == 1
        % plot(locMeanCorr(np),PCTloc0(np),'.','Color',colors.PV); % PV
        plot(locMeanCorr1(np),PCTloc1(np),'.','Color',[0.7 0 0]); % PV
        locCorrPV = [locCorrPV; locMeanCorr1(np)];
        locAccPV = [locAccPV; PCTloc1(np)];
    elseif locInhType(np) == 2
        % plot(locMeanCorr(np),PCTloc0(np),'.','Color',colors.SOM); % SOM
        plot(locMeanCorr1(np),PCTloc1(np),'.','Color',[0.7 0 0]); % SOM
        locCorrSOM = [locCorrSOM; locMeanCorr1(np)]; 
        locAccSOM = [locAccSOM; PCTloc1(np)];
    elseif locInhType(np) == 0
        % only EXC neurons present
        plot(locMeanCorr1(np),PCTloc1(np),'.','Color',[0.7 0 0]); 
    end
end

% SAVE AVERAGE ACCURACIES
accMat(1,3) = mean(PCTloc1,'omitnan');
allAcc{3} = PCTloc1;
locCI = mkCI(PCTloc1);
accMat(2,3) = accMat(1,3) - locCI(1);
accMatAA(1,3) = mean(PCTloc1AA,'omitnan');
allAccAA{3} = PCTloc1AA;
locCI = mkCI(PCTloc1AA);
accMatAA(2,3) = accMatAA(1,3) - locCI(1);

% FIT ALL POINTS
mdl = fitlm(locMeanCorr1,PCTloc1);
xGrid = linspace(min(locMeanCorr1(~isnan(PCTloc1))), max(locMeanCorr1(~isnan(PCTloc1))), 100)'; 
[y,ci] = predict(mdl,xGrid);
plot(xGrid, y, 'r-', xGrid, ci(:,1), 'r:', xGrid, ci(:,2), 'r:', 'LineWidth', 1);
disp("EXC:EXC base window fit");
disp(strcat("p-val = ",num2str(mdl.ModelFitVsNullModel.Pvalue,3)));
disp(strcat("r-sq. = ",num2str(mdl.Rsquared.Ordinary,3)));
disp(strcat("slope = ",num2str(mdl.Coefficients.Estimate(2),3)));



% EXC:EXC LATE WINDOW: (3-4sec post-stim) 
PCTloc = squeeze(mean(PCTallEXC_late(:,~trainedSess,1),1))' .* 100;
PCTloc(PCTloc == 0) = NaN;

% split pairwise correlations up by cell type (PV|SOM|EXC)
allPairTypes = cell(size(allCorrValsLate));
for ns = 1 : size(allPairNumsLate,1)
    allPairTypes{ns} = zeros(length(allPairNumsLate{ns}),1);
    for np = 1 : length(allPairNumsLate{ns})
        locPairID = Tc.identity(allPairNumsLate{ns}(np,:));
        if all(strcmp(locPairID,'EXC')) % EXC & EXC [1]
            allPairTypes{ns}(np) = 1;
        elseif any(strcmp(locPairID,'EXC')) && any(strcmp(locPairID,'PV')) % EXC & PV [2]
            allPairTypes{ns}(np) = 2;
        elseif any(strcmp(locPairID,'EXC')) && any(strcmp(locPairID,'SOM')) % EXC & SOM [3]
            allPairTypes{ns}(np) = 3;
        elseif all(strcmp(locPairID,'PV')) % PV & PV [4]
            allPairTypes{ns}(np) = 4;
        elseif all(strcmp(locPairID,'SOM')) % SOM & SOM [5]
            allPairTypes{ns}(np) = 5;
        end
    end
end

% EXTRACT CORRELATION VALUES FROM WITHIN SESSION FREQUENCY-SELECTIVE PAIRS 
npt = 1;
locPairVals = NaN(size(allPairNumsLate,1),1);
locMeanCorr2 = NaN(size(allPairNumsLate,1),1);
locInhType = zeros(size(allPairNumsLate,1),1);
for ns = 1 : size(allPairNumsLate,1) % loop sessions
    if multiAmp == true
        locIdx = all(is.allMod(allPairNumsLate{ns}),2); % neuron pairs where both are selective
    elseif multiAmp == false
        locIdx = all(is.selective(allPairNumsLate{ns}),2); % neuron pairs where both are selective
    end
    locVals = allCorrValsLate{ns}(allPairTypes{ns} == npt & locIdx);
    locPairVals(ns) = mean(locVals);
    locMeanCorr2(ns) = mean(locVals); % corr values from selective pairs
    if any(allPairTypes{ns}==2) % PV present
        locInhType(ns) = 1;
    elseif any(allPairTypes{ns}==3) % SOM present
        locInhType(ns) = 2;
    end
end
locMeanCorr2 = locMeanCorr2(~trainedSess);

% APPLY THRESHOLD TO WINTIN-SESSION AVG NEURON CORRELATION VALUES
if multiAmp == true
    abvThresh = locMeanCorr2 > 0.3;
elseif multiAmp == false
    abvThresh = locMeanCorr2 > 1;
end
locMeanCorr2(abvThresh) = [];
PCTloc2 = PCTloc(~abvThresh);
locInhType(abvThresh) = [];

% PLOTTING
locCorrPV = []; locAccPV = [];
locCorrSOM = []; locAccSOM = [];
for np = 1 : length(locMeanCorr2)
    if locInhType(np) == 1
        % plot(locMeanCorr(np),PCTloc0(np),'.','Color',colors.PV); % PV
        plot(locMeanCorr2(np),PCTloc2(np),'.','Color',[0 0 0.7]); % PV
        locCorrPV = [locCorrPV; locMeanCorr2(np)]; 
        locAccPV = [locAccPV; PCTloc2(np)];
    elseif locInhType(np) == 2
        % plot(locMeanCorr(np),PCTloc0(np),'.','Color',colors.SOM); % SOM
        plot(locMeanCorr2(np),PCTloc2(np),'.','Color',[0 0 0.7]); % SOM
        locCorrSOM = [locCorrSOM; locMeanCorr2(np)]; 
        locAccSOM = [locAccSOM; PCTloc2(np)];
    elseif locInhType(np) == 0
        % only EXC nuerons present in sess
        plot(locMeanCorr2(np),PCTloc2(np),'.','Color',[0 0 0.7]); 
    end
end

% FIT ALL POINTS
mdl = fitlm(locMeanCorr2,PCTloc2);
xGrid = linspace(min(locMeanCorr2(~isnan(PCTloc2))), max(locMeanCorr2(~isnan(PCTloc2))), 100)'; 
[y,ci] = predict(mdl,xGrid);
plot(xGrid, y, 'b-', xGrid, ci(:,1), 'b:', xGrid, ci(:,2), 'b:', 'LineWidth', 1);
disp(strcat("EXC:EXC late window fit p-val = ",num2str(mdl.ModelFitVsNullModel.Pvalue,3)));
disp(strcat("EXC:EXC late window fit r-sq. = ",num2str(mdl.Rsquared.Ordinary,3)));
disp(strcat("slope = ",num2str(mdl.Coefficients.Estimate(2),3)));

% FINAL FIGURE ADJUSTMENTS
ylabel('Frequency decoding accuracy %','FontName','Arial');
xlabel('Mean neuron corr value','FontName','Arial');
title({"Neuron correlation by frequency decoding accuracy",pairNames{npt}},'FontName','Arial');
set(gca,'XLim',[.1 .8],'YLim',[10 80]);
hold off;

% ADDITIONAL SUB PLOT OF CORRELATION SCATTER
figure('Color',[1 1 1]);
plot(locMeanCorr1,locMeanCorr2,'k.');
set(gca,'XLim',[0 .8],'YLim',[0 .8],'box','off');



% EXC:PV and EXC:SOM
for npt = 2 : 3 % num pair types (EXC:EXC, EXC:PV, EXC:SOM, PV:PV, SOM:SOM)
    figure('Color',[1 1 1]); hold on;
    
    % BASELINE WINDOW (0-1 sec post-stim)
    PCTloc = squeeze(mean(PCTall_baseline(:,~trainedSess,1),1))' .* 100;
    PCTloc(PCTloc == 0) = NaN;
    PCTlocAA = squeeze(mean(PCTall_baseline_allAvail(:,~trainedSess,1),1))' .* 100;
    PCTlocAA(PCTlocAA == 0) = NaN;

    % PREALLOCATE
    locPairVals = NaN(size(allPairNumsBase,1),1);
    locMeanCorr1 = NaN(size(allPairNumsBase,1),1);
    for ns = 1 : size(allPairNumsBase,1) % loop sessions
        if multiAmp == true
            locIdx = all(is.allMod(allPairNumsBase{ns}),2); % neuron pairs where both are selective
        elseif multiAmp == false
            locIdx = all(is.selective(allPairNumsBase{ns}),2); % neuron pairs where both are selective
        end
        locVals = allCorrValsBase{ns}(allPairTypes{ns} == npt & locIdx);
        locPairVals(ns) = mean(locVals);
        locMeanCorr1(ns) = mean(locVals); % corr values from selective pairs
    end
    locMeanCorr1 = locMeanCorr1(~trainedSess);
    PCTloc(isnan(locMeanCorr1)) = NaN;
    PCTlocAA(isnan(locMeanCorr1)) = NaN;

    % SAVE OUT AVERAGE ACCURACIES
    accMat(1,npt-1) = mean(PCTloc,'omitnan'); % FREQ.SEL.NEURONS
    allAcc{npt-1} = PCTloc;
    locCI = mkCI(PCTloc);
    accMat(2,npt-1) = accMat(1,npt-1) - locCI(1);
    accMatAA(1,npt-1) = mean(PCTlocAA,'omitnan'); % ALL.AVAIL.NEURONS
    allAccAA{npt-1} = PCTlocAA;
    locCI = mkCI(PCTlocAA);
    accMatAA(2,npt-1) = accMatAA(1,npt-1) - locCI(1);
    plot(locMeanCorr1,PCTloc,'.','Color',[0.7 0 0]);

    % LINEAR MODEL FIT
    mdl = fitlm(locMeanCorr1,PCTloc);
    xGrid = linspace(min(locMeanCorr1(~isnan(PCTloc1))), max(locMeanCorr1(~isnan(PCTloc1))), 100)'; 
    [y,ci] = predict(mdl,xGrid);
    plot(xGrid, y, 'r-', xGrid, ci(:,1), 'r:', xGrid, ci(:,2), 'r:', 'LineWidth', 1);
    disp(strcat(pairNames{npt},"base window fit p-val = ",num2str(mdl.ModelFitVsNullModel.Pvalue,3)));
    disp(strcat(pairNames{npt},"base window fit r-sq. = ",num2str(mdl.Rsquared.Ordinary,3)));
    disp(strcat("slope = ",num2str(mdl.Coefficients.Estimate(2),3)));

    % LATE WINDOW (3-4 sec post-stim)
    PCTloc = squeeze(mean(PCTall_late(:,~trainedSess,1),1))' .* 100;
    PCTloc(PCTloc == 0) = NaN;
    locPairVals = NaN(size(allPairNumsLate,1),1);
    locMeanCorr2 = NaN(size(allPairNumsLate,1),1);
    for ns = 1 : size(allPairNumsLate,1) % loop sessions
        if multiAmp == true
            locIdx = all(is.allMod(allPairNumsLate{ns}),2); % neuron pairs where both are selective
        elseif multiAmp == false
            locIdx = all(is.selective(allPairNumsLate{ns}),2); % neuron pairs where both are selective
        end
        locVals = allCorrValsLate{ns}(allPairTypes{ns} == npt & locIdx);
        locPairVals(ns) = mean(locVals);
        locMeanCorr2(ns) = mean(locVals); % corr values from selective pairs
    end
    locMeanCorr2 = locMeanCorr2(~trainedSess);
    plot(locMeanCorr2,PCTloc,'.','Color',[0 0 0.7]);

    % LINEAR MODEL FIT
    mdl = fitlm(locMeanCorr2,PCTloc);
    xGrid = linspace(min(locMeanCorr2(~isnan(PCTloc2))), max(locMeanCorr2(~isnan(PCTloc2))), 100)'; 
    [y,ci] = predict(mdl,xGrid);
    plot(xGrid, y, 'b-', xGrid, ci(:,1), 'b:', xGrid, ci(:,2), 'b:', 'LineWidth', 1);
    disp(strcat(pairNames{npt}," late window fit p-val = ",num2str(mdl.ModelFitVsNullModel.Pvalue,3)));
    disp(strcat(pairNames{npt}," late window fit r-sq. = ",num2str(mdl.Rsquared.Ordinary,3)));
    disp(strcat("slope = ",num2str(mdl.Coefficients.Estimate(2),3)));

    % FIGURE ADJUSTMENTS 
    set(gca,'XLim',[.1 .8],'YLim',[10 80]);
    ylabel('Frequency decoding accuracy %','FontName','Arial');
    xlabel('Mean neuron corr value','FontName','Arial');
    title({"Neuron correlation by frequency decoding accuracy",pairNames{npt}},'FontName','Arial');
    if npt == 3
        legend({"dF/F window (0-1sec)","linear fit","95% CI","","dF/F window (3-4sec)","linear fit","95% CI",""},...
            'Box','off','Location','NorthWest');
    end
    hold off;

    % ADDITIONAL SUBPLOT | CORR SCATTER
    figure('Color',[1 1 1]);
    plot(locMeanCorr1,locMeanCorr2,'k.');
    set(gca,'XLim',[0 .8],'YLim',[0 .8],'box','off');

end




%%  Additional plost: average accuracy across 3 conditions (EXC:EXC, EXC:PV, EXC:SOM)
% FREQUENCY-SELECTIVE UNITS ONLY
figure('Color',[1 1 1]); hold on;
for nc = 1 : length(cellType)
    bar(nc,accMat(1,nc),'FaceColor',colors.(cellType(nc)));
    errorbar(nc,accMat(1,nc),accMat(2,nc),'vertical','Color',[0 0 0]);
end
set(gca,'XLim',[0.2 3.8],'XTick',[2 6],'XTickLabel',{'Frequency selective','all available'});
ylabel('Mean frequency decoding accuracy %','FontName','Arial');
xlabel('Neuron source(s)','FontName','Arial');
legend({'PV+EXC','','SOM+EXC','','EXC alone',''},'Box','off');
title("Average within-session decoding accuracy",'FontName','Arial');

% ALL AVAILABLE UNITS
figure('Color',[1 1 1]); hold on;

% PV+EXC BAR
yVals1 = allAcc{1}(~isnan(allAcc{1}));
yMean = mean(yVals1);
yCI = mkCI(yVals1);
bar(1,yMean,'FaceColor',colors.PV);
errorbar(1,yMean,yMean-yCI(1),'vertical','Color',[0 0 0]);
% EXC-ALONE BAR (NO PV)
yVals2 = allAcc{3}(~isnan(allAcc{1}));
yMean = mean(yVals2,'omitnan');
yCI = mkCI(yVals2);
bar(2,yMean,'FaceColor',colors.EXC);
errorbar(2,yMean,yMean-yCI(1),'vertical','Color',[0 0 0]);
% 2-sided Kolmogorov-Smirnov test [PV+EXC : EXC-alone]
[~,p,~] = kstest2(yVals1,yVals2);
disp(p);
plot([1 2],[53 53],'k-');
if p < 0.05
    text(1.5,53,'*','HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontName','Arial','FontSize',12);
else
    text(1.5,53,'n.s.','HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontName','Arial','FontSize',12);
end
% SOM+EXC BAR
yVals3 = allAcc{2}(~isnan(allAcc{2}));
yMean = mean(yVals3);
yCI = mkCI(yVals3);
bar(3,yMean,'FaceColor',colors.SOM);
errorbar(3,yMean,yMean-yCI(1),'vertical','Color',[0 0 0]);
% EXC-ALONE BAR (NO SOM)
yVals4 = allAcc{3}(~isnan(allAcc{2}));
yMean = mean(yVals4,'omitnan');
yCI = mkCI(yVals4);
bar(4,yMean,'FaceColor',colors.EXC);
errorbar(4,yMean,yMean-yCI(1),'vertical','Color',[0 0 0]);
% 2-sided Kolmogorov-Smirnov test [PV+EXC : EXC-alone]
[~,p,~] = kstest2(yVals3,yVals4);
disp(p);
plot([3 4],[53 53],'k-');
if p < 0.05
    text(3.5,53,'*','HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontName','Arial','FontSize',12);
else
    text(3.5,53,'n.s.','HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontName','Arial','FontSize',12);
end

% FIGURE ADJUSTMENTS
set(gca,'XLim',[0.2 4.8],'XTick',1:4,'YLim',[32.5 57.5],'XTickLabel',{'PV+EXC','EXC alone','SOM+EXC','EXC alone'},...
    'XTickLabelRotation',-20);
ylabel('Mean frequency decoding accuracy %','FontName','Arial');
xlabel('Neuron source(s)','FontName','Arial');
title("Average within-session decoding accuracy",'FontName','Arial');



%% PV:EXC vs SOM:EXC decoding as histograms

% PV+EXC BAR
yVals1 = allAcc{1}(~isnan(allAcc{1}));

% SOM+EXC BAR
yVals3 = allAcc{2}(~isnan(allAcc{2}));

% ALL AVAILABLE UNITS
figure('Color',[1 1 1]); hold on;
histogram(yVals1,'FaceAlpha',0.8,'FaceColor',colors.PV);
histogram(yVals3,'FaceAlpha',0.8,'FaceColor',colors.SOM);








