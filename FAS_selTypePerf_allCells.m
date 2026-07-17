% FAS_selTypePerf_allCells.m
%
% Selectivity Type Performance from all cells
%
% Dataset2, selective draws vs. random draws

%% Preload model results

%  cd('C:\Users\skich\Desktop');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% % [5 amps, 200 iterations, 4 sel types, 15:5:N draws, 3 cell types]
% tic; load("C:\Users\skich\Desktop\WORK\modelResults_PCTall_MDLall.mat"); toc;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% % 5 amps, 200 iterations, 1 sel type (all freq), 15:5:N draws, 3 cell types]
% tic; load("C:\Users\skich\Desktop\WORK\modelResults_PCTall2_MDLall2.mat"); toc;
% tic; load("C:\Users\skich\Desktop\WORK\randDraws.mat"); toc;
% load("C:\Users\skich\Desktop\WORK\modelAccuracy_dataset2_allFreqSel.mat");
% PCTall = cat(3,PCTall1,PCTall2);
% clear PCTall1 PCTall2

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% % [top amp, 200 iterations, 4 sel types, 15:5:990 draws, 3EXC]
% tic; load("C:\Users\skich\Desktop\modelResults_PCTall3_MDLall3.mat"); toc;

% DATA: FreqSel VS All neurons VS random shuffle (2025-11-19)
cd('C:\Users\skich\Desktop\WORK');
load('modelAccuracy_dataset2_allFreqSel.mat');
% [1] Amplitude
% [2] Repetitions
% [3] Source (FreqSel, AllNeurons, RandomShuffleData)
% [4] Bin number (15:5:300)
% [3] CellType

% DATA: 4 SEL TYPES (2025-11-19)
cd('C:\Users\skich\Desktop\WORK');
load('modelAccuracy_dataset2_allSelTypes.mat');


%% perform population decoding

nIter = 200; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
targetAmp = 5; % 0.33,0.76,1.6,3.5,7.7
usePCA = true; % use PCA step prior to model training
kStart = 1; % first PC to include
kComp = 15; % total number of PCs to use in model training / testing
useFreqEffNeu = true; % use neurons with significant anova freq effects
useAllNeu = false; % use all available neurons for modeling
includeErrBar = true; %include error bars in plots
cellType = ["PV","SOM","EXC"]; % order for plotting
% cellType = "EXC"; % order for plotting
% selType = ["FreqMod","AmpMod","DualMod","InterMod"]; 
is.allMod = (is.FreqMod | is.DualMod | is.InterMod) & ~is.trained;
selType = "allMod";
% is.all = true(nCellTotal,1);
% selType = "all";
neuStepSize = 5; % num neurons to add to model each step
numCellDraws = 1; % over-estimation to poulate arrays

if length(selType) <= 2
    maxDraw = 172; % 2 SelTypes: FreqVsRandom (300), 
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
PCTall = zeros(nAmp,nIter,length(selType),numCellDraws,length(cellType)); % (na,nj,nt,nd,nc)
% Models (all)
% MDLall = cell(nIter,length(selType),length(nCellDraws),length(cellType)); % (nj,nt,nd,nc)
MDLall = cell(nAmp,nIter,length(selType),numCellDraws,length(cellType)); % (na,nj,nt,nd,nc)

predFreqAll = cell(nFreq,nIter,length(selType),numCellDraws,length(cellType));
testFreqAll = cell(nFreq,nIter,length(selType),numCellDraws,length(cellType));

% % Create output structure for storing PCA information
% PC = struct();
% for nc = 1 : length(cellType)
%     PC.(cellType(nc)) = struct();
% end


for na = 1 : nAmp %numAmps
targetAmp = na;
% na = targetAmp; %only use one amplitude
disp(strcat("Amp: ",num2str(ampValAll(na))));

for nc = 1:length(cellType) %nt = 1 : length(selType)
disp(cellType(nc));
% disp(selType(nt));

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
    if useFreqEffNeu == false
        maxNeuAvail = 120;
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
    % predFreqAll = zeros(nIter,nCellType,nFreq*testTrial);
    % testFreqAll = zeros(nIter,nCellType,nFreq*testTrial);
    
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
                [coeff,score,latent] = pca(trainResp);
                % PC.(cellType(nc)).coeff{ni} = coeff;
                % PC.(cellType(nc)).score{ni} = score;
                % PC.(cellType(nc)).latent{ni} = latent;

                trainResp_red = score(:,kStart:kComp);
                mdl = fitcdiscr(trainResp_red,trainFreq,'DiscrimType','linear'); %train model
                % mdl = fitcdiscr(trainResp_red,trainFreq(randperm(length(trainFreq))),'DiscrimType','linear'); %train model (RAND MIXED DATA)
                predFreq = predict(mdl,testResp*coeff(:,kStart:kComp));
            else
                % WITHOUT PCA
                mdl = fitcdiscr(trainResp,trainFreq,'DiscrimType','linear');
                predFreq = predict(mdl,testResp);
            end
    
            % collect output
            PCTall(na,nj,nt,nd,nc) = mean(predFreq == testFreq); % accuracy
            MDLall{na,nj,nt,nd,nc} = mdl; % model parameters
            predFreqAll{na,nj,nt,nd,nc} = predFreq;
            testFreqAll{na,nj,nt,nd,nc} = testFreq;

        end 
    end
    fprintf(newline)
end

end

end % na (numAmps)


%% PLOT: PERFORMANCE BY CELL TYPE
% TOP AMP, ALL FREQ NEURONS POOLED

% KEY: 
% [nj] Model Iterations     ... [scalar, 200]
% [nt] Selectivity Types    ... ["FreqMod","AmpMod","DualMod","InterMod"]
% [nd] Number of Cell Draws ... [various lengths, 5-8 total]
% [nc] Cell Types           ... ["EXC","PV","SOM"]

targetAmp = 5;
neuStepSize = 5;
selType = "allMod";
%cellType = ["PV","SOM","EXC"];
maxDraw = 300;
figure('Color',[1 1 1],'Theme','Light'); hold on;

% fake data for correct legend order
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');
plot([-1 -2],[-1 -1],'k-');  % freq sel
plot([-1 -2],[-1 -1],'k--'); % all neurons
% plot([-1 -2],[-1 -1],'k:'); % rand shuffle

% PLOT DATA FROM PCTall2 (FREQ. SEL NEURONS BINNED)
for nt = 1 : 2 %size(PCTall,3)
    for nc = 1 : length(cellType)
        switch nt
            case 1
                if strcmp(cellType(nc),'EXC')
                    maxNeuAvail = maxDraw;
                else % PV | SOM
                    maxNeuAvail = sum(is.(cellType(nc)) & is.(selType(nt)));
                end
            case 2
                maxNeuAvail = 300;
        end
        nCellDraws = 15:neuStepSize:maxNeuAvail;
        target = cellType(nc);
        nBins = 1:length(nCellDraws);
        yMean = mean(squeeze(PCTall(targetAmp,:,nt,nBins,nc)),1).*100;
        % ySems = sem(squeeze(PCTall2(targetAmp,:,nt,nBins,nc)),1).*100;
        yCIs = zeros(length(nBins),2); % confidence intervals (95%)
        for nd = 1 : length(nBins)
            yCIs(nd,:) = mkCI(squeeze(PCTall(targetAmp,:,nt,nd,nc))).*100;
        end
        xLooped = [nCellDraws flip(nCellDraws)]; 
        % yLooped = [yMean+ySems flip(yMean-ySems)]; 
        yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
        fill(xLooped,yLooped,colors.(target),'FaceAlpha',0.3,'EdgeColor','none');
        switch nt 
            case 1       
                plot(nCellDraws,yMean,'-','LineWidth',1,'Color',colors.(target));
            case 2
                plot(nCellDraws,yMean,'--','LineWidth',1,'Color',colors.(target));
            case 3
                plot(nCellDraws,yMean,':','LineWidth',1,'Color',colors.(target));
        end
    end
end

ylabel('model accuracy (%)','FontName','Arial');
xlabel('Neurons used','FontName','Arial');
title('Frequency decoding performance by cell type','FontName','Arial');
subTitleTxt = strcat("From trials at amplitude: ",num2str(ampValAll(targetAmp)),"\mum");
subtitle(subTitleTxt,'FontName','Arial');
% set(gca,'XLim',[13 300],'YLim',[18 82],'YTick',10:10:100,'XTick',15:30:300);
set(gca,'XLim',[13 300],'YLim',[20 82],'YTick',10:10:100,'XTick',15:30:300);
legend({"PV","SOM","EXC","Freq. selective","All available"},'Location','Northwest',...
    'box','off','NumColumns',2);


%% FREQUENCY DECODING PERFORMANCE ACROSS AMPLITUDES

% cd('C:\Users\skich\Desktop\WORK');
% load('modelAccuracy_dataset2_allFreqSel.mat');

targetDrawBin = 22; % bin 30 is 160 neurons (15:5:300)
xVals = ampValAll'; 
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

% plot([15 maxDrawEXC],[20 20],'r--'); % chance performance level
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Amplitude (\mum)','FontName','Arial');
title('Frequency decoding performance across amplitudes','FontName','Arial');
subtitle('Fixed model size: 160 frequency selective neurons','FontName','Arial');
set(gca,'XLim',[0.32 7.9],'YLim',[35 75],'YTick',10:5:100,'XTick',xVals,'XTickLabelRotation',45,...
    'XScale','log');
% legend({"Freq. Mod","All neurons","EXC","PV","SOM"},'Location','Northwest');
% legend({"PV","SOM","EXC","Freq. selective","All available"},'Location','Northwest',...
%     'box','off','NumColumns',2);
legend({"PV","SOM","EXC"},'Location','Northwest',...
    'box','off','NumColumns',1);



%% PLOT PC TUNING CURVES AS HEAT MAPS
targetDrawBin = 18;
targetRep = 1;
targetAmp = 5;
numPCs = 15;
RGB = mkcolors([1 0 0],[0 0 1],numPCs);
figure('Theme','light','color',[1 1 1]);
tiledlayout(1,length(cellType));
for nc = 1 : length(cellType)
    mdl = MDLall2{targetAmp,targetRep,1,targetDrawBin,nc}; 
    PCresps = zeros(12,5,numPCs);
    nexttile; hold on;
    for np = 1 : numPCs
        for nf = 1 : 5
            PCresps(:,nf,np) = mdl.X(mdl.Y == freqVal(nf),np);
            %ySem = sem(PCresps(:,nf,np));
            yMean = mean(PCresps(:,nf,np));
            %plot([nf nf],[yMean-ySem yMean+ySem],'-','Color',RGB(np,:)); %sem lines
        end
        %plot(1:5, mean(squeeze(PCresps(:,:,np)),1),'-','Color',RGB(np,:));%mean line
    end
    imagesc(flip(squeeze(mean(PCresps,1))'));
    subtitle(cellType(nc),'FontName','Arial');
    set(gca,'YTick',[1 6 11 15],'YTickLabel',{"15","10","5","1"},'XTick',1:5,...
        'XTickLabel',freqVal);
    if nc == 1 
        ylabel('PC Number','FontName','Arial');
    end
    if nc == 2
        title("Frequency tuning of PCs",'FontName','Arial');
        xlabel('frequency (Hz)','FontName','Arial');
    end
    hold off
end

%% PLOT PC TUNING CURVES FROM MODELS
targetCellType = 3; %(PV, SOM, EXC)
targetDrawBin = 1;
targetAmp = 5;
RGB = mkcolors([1 0 0],[0 0 1],10);
figure('Theme','light','color',[1 1 1]);
tiledlayout(2,2);
for nr = 1 : 4 %num reps
mdl = MDLall{targetAmp,nr,1,targetDrawBin,targetCellType}; % nIter (200), nCellTypes (3), nCellDraws (4)
PCresps = zeros(12,5,15);
nexttile(nr); hold on;
for np = 1 : 7
    for nf = 1 : 5
        PCresps(:,nf,np) = mdl.X(mdl.Y == freqVal(nf),np);
        ySem = sem(PCresps(:,nf,np));
        yMean = mean(PCresps(:,nf,np));
        plot([nf nf],[yMean-ySem yMean+ySem],'-','Color',RGB(np,:)); %sem lines
    end
    plot(1:5, mean(squeeze(PCresps(:,:,np)),1),'-','Color',RGB(np,:));%mean line
end
set(gca,'XLim',[0.5 5.5],'XTick',1:5,'XTickLabel',freqVal);
ylabel('z-scored "response"','FontName','Arial');
xlabel('frequency (Hz)','FontName','Arial');
hold off;
end


%% PLOT: PERFORMANCE BY CELL TYPE AND SEL TYPE
% TOP AMP, EACH SELECTIVITY TYPE SEPARATE

% KEY: 
% [nj] Model Iterations     ... [scalar, 200]
% [nt] Selectivity Types    ... ["FreqMod","AmpMod","DualMod","InterMod"]
% [nd] Number of Cell Draws ... [various lengths, 5-8 total]
% [nc] Cell Types           ... ["EXC","PV","SOM"]

targetAmp = 1;
maxDraw = 120;
neuStepSize = 5;
selType = ["FreqMod","AmpMod","DualMod","InterMod"]; 
figure('Color',[1 1 1],'Theme','Light'); hold on;
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
        yMean = mean(squeeze(PCTall(targetAmp,:,nt,nBins,nc)),1).*100;
        ySems = sem(squeeze(PCTall(targetAmp,:,nt,nBins,nc)),1).*100;
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
xlabel('Neurons used','FontName','Arial');
title('Frequency decoding accuracy by selectivity type','FontName','Arial');
set(gca,'XLim',[15 maxDraw],'YLim',[20 80],'YTick',20:10:80,'XTick',15:15:120);
plot([15 maxDraw],[20 20],'r--');

% ADDITIONAL FIG FOR LEGEND
figure('theme','light','color',[1 1 1]); hold on;
plot([1 2],[1 1],'k-'); plot([1 2],[1 1],'k-.'); % Inter & Dual
plot([1 2],[1 1],'k--'); plot([1 2],[1 1],'k:'); % Freq & Amp
fill([1 2 2 1],[1.1 1.1 .9 .9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
fill([1 2 2 1],[1.1 1.1 .9 .9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
fill([1 2 2 1],[1.1 1.1 .9 .9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');
plot([1 2],[1 1],'r--'); % chance performance level
set(gca,'XLim',[1 5],'YLim',[1 5]);
legend({"Interaction","Both mod.","Freq. mod.","Amp. mod.","PV","SOM","EXC","chance"},'NumColumns',2);




%% PLOT: PERFORMANCE BY CELL TYPE AND SEL TYPE
% TOP AMP, EACH SELECTIVITY TYPE SEPARATE
% UPDATED 2026/01/26

% KEY: 
% [nj] Model Iterations     ... [scalar, 200]
% [nt] Selectivity Types    ... ["FreqMod","AmpMod","DualMod","InterMod"]
% [nd] Number of Cell Draws ... [various lengths, 5-8 total]
% [nc] Cell Types           ... ["EXC","PV","SOM"]

cellType = ["PV","SOM","EXC"];
cellType = "EXC";
targetAmp = 5;
maxDraw = 120;
neuStepSize = 5;
selType = ["FreqMod","AmpMod","DualMod","InterMod"]; 
figure('Color',[1 1 1],'Theme','Light'); hold on;
for nt = 1 : length(selType)
    for nc = 1 : length(cellType)
        if strcmp(cellType(nc),'EXC')
            maxNeuAvail = maxDraw;
        else % PV | SOM
            maxNeuAvail = sum(is.(cellType(nc)) & is.(selType(nt)));
            idxNeuBin = squeeze(PCTall(1,1,1,:,nc))~=0;
        end
        nCellDraws = 15:neuStepSize:maxNeuAvail;
        target = cellType(nc);
        nBins = 1:length(nCellDraws);
        yMean = mean(squeeze(PCTall(targetAmp,:,nt,nBins,nc)),1).*100;
        ySems = sem(squeeze(PCTall(targetAmp,:,nt,nBins,nc)),1).*100;
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
xlabel('Neurons used','FontName','Arial');
title('Frequency decoding accuracy by selectivity type','FontName','Arial');
set(gca,'XLim',[15 maxDraw],'YLim',[20 84],'YTick',20:10:100,'XTick',15:15:120);
plot([15 maxDraw],[20 20],'r--');

% ADDITIONAL FIG FOR LEGEND
figure('theme','light','color',[1 1 1]); hold on;
plot([1 2],[1 1],'k-'); plot([1 2],[1 1],'k-.'); % Inter & Dual
plot([1 2],[1 1],'k--'); plot([1 2],[1 1],'k:'); % Freq & Amp
fill([1 2 2 1],[1.1 1.1 .9 .9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
fill([1 2 2 1],[1.1 1.1 .9 .9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
fill([1 2 2 1],[1.1 1.1 .9 .9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');
plot([1 2],[1 1],'r--'); % chance performance level
set(gca,'XLim',[1 5],'YLim',[1 5]);
legend({"Interaction","Both mod.","Freq. mod.","Amp. mod.","PV","SOM","EXC","chance"},'NumColumns',2);



%% ACCURACY AT N=40, EACH OF 4 CELL TYPES
targetAmp = 5; % 7.7 um
targetBin = 6; % 40 cell draws
xSpace = [1,2,3;6,7,8;11,12,13;16,17,18];
ntOrder = [2 1 3 4];
ntNames = {"A","F","F+A","F*A"};
figure('theme','light','color',[1 1 1]); hold on;
plot(-1,-1,'.','MarkerSize',8,'Color',colors.PV);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.SOM);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.EXC);
for nt = 1 : length(selType)
    for nc = 1 : length(cellType)
        yVals = PCTall(targetAmp,:,ntOrder(nt),targetBin,nc);
        yMean = mean(yVals)*100; % mean y value fonverted to percentage
        yCIs = (mkCI(yVals) - mean(yVals)).*100; % 95% confidence intervals centered at mean
        plot(xSpace(nt,nc),yMean,'.','MarkerSize',8,'Color',colors.(cellType(nc)));
        errorbar(xSpace(nt,nc),yMean,yCIs(1),yCIs(2),'Color',colors.(cellType(nc)));
    end
end

ylabel('model accuracy (%)','FontName','Arial');
xlabel('Selectivity Type','FontName','Arial');
title('Frequency decoding performance of each selectivity type','FontName','Arial');
subtitle('Amplitude: 0.33\mum, model size: 40 neurons ea.','FontName','Arial');
set(gca,'YLim',[15 65],'XLim',[0 18],'XTick',xSpace(:,2),'XTickLabel',ntNames,'FontName','Arial');
legend({"PV","SOM","EXC"},'Location','Northwest','FontName','Arial','box','off');



%% ACCURACY AT N=40, EACH OF 4 CELL TYPES, EACH INDV  AMP TOGETHER

targetBin = 6; % 40 cell draws
nAmps = size(PCTall,1);
xSpace = [1:5; 8:12; 15:19; 22:26];

ntOrder = [2 1 3 4];
ntNames = {"Amplitude","Frequency","Dual Effect","Interaction"};
figure('theme','light','color',[1 1 1]); hold on;
plot(-1,-1,'.','MarkerSize',8,'Color',colors.PV);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.SOM);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.EXC);


for nc = 1 : length(cellType)
    yMeans = zeros(4,nAmps);
    for nt = 1 : length(ntNames)
        yVals=PCTall(:,:,ntOrder(nt),targetBin,nc);
        yMeans(nt,:) = mean(yVals,2).*100; % mean y value converted to percentage
        yCIs = zeros(size(yVals,1),2);
        for na = 1 : nAmps
            yCIs(na,:) = (mkCI(yVals(na,:)) - mean(yVals(na,:))).*100; % 95% confidence intervals centered at mean
        end
        plot(xSpace(nt,:),yMeans(nt,:),'-','MarkerSize',8,'Color',colors.(cellType(nc)));
        plot(xSpace(nt,:),yMeans(nt,:),'.','MarkerSize',8,'Color',colors.(cellType(nc)));
    end
end

ylabel('model accuracy (%)','FontName','Arial');
xlabel('Selectivity Type','FontName','Arial');
title('Frequency decoding accuracy by selectivity type','FontName','Arial');
subtitle('40 neurons per condition','FontName','Arial');
set(gca,'YLim',[17 70],'XLim',[0 27],'XTick',xSpace(:,3),'XTickLabel',ntNames,'FontName','Arial');
legend({"PV","SOM","EXC"},'Location','Northwest','FontName','Arial',...
    'box','off');


%% ACCURACY AT N=40, EACH OF 4 CELL TYPES, COMBINE ACROSS AMPS
targetBin = 6; % 40 cell draws
nAmps = size(PCTall,1);
xSpace = [1,2,3;7,8,9;13,14,15;19,20,21];
xGap = 0.5;
ntOrder = [2 1 3 4];
ntNames = {"Amplitude","Frequency","Dual Effect","Interaction"};
figure('theme','light','color',[1 1 1]); hold on;
plot(-1,-1,'.','MarkerSize',8,'Color',colors.PV);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.SOM);
plot(-1,-1,'.','MarkerSize',8,'Color',colors.EXC);


for nc = 1 : length(cellType)
    yMeans = zeros(4,nAmps);
    for nt = 1 : length(ntNames)
        %for na = 1 : length(targetAmp)
            yVals=PCTall(:,:,ntOrder(nt),targetBin,nc);
            yMeans(nt,:) = mean(yVals,2).*100; % mean y value converted to percentage
            yMean = mean(yVals,'all')*100;
            yCI = mkCI(yVals(:)).*100;
            yCIs = zeros(size(yVals,1),2);
            for na = 1 : nAmps
                yCIs(na,:) = (mkCI(yVals(na,:)) - mean(yVals(na,:))).*100; % 95% confidence intervals centered at mean
            end
            plot(xSpace(nt,nc),yMean,'.','MarkerSize',8,'Color',colors.(cellType(nc)));
            errorbar(xSpace(nt,nc),yMean,yCI(1)-yMean,yCI(2)-yMean,'Color',colors.(cellType(nc)));
        %end
    end
end


ylabel('model accuracy (%)','FontName','Arial');
xlabel('Selectivity Type','FontName','Arial');
title('Frequency decoding accuracy by selectivity type','FontName','Arial');
subtitle('40 neurons per condition','FontName','Arial');
set(gca,'YLim',[17 53],'XLim',[0 22],'XTick',xSpace(:,2)+(xGap/2),'XTickLabel',ntNames,'FontName','Arial');
legend({"PV","SOM","EXC"},'Location','Northwest','FontName','Arial',...
    'box','off');

%% ACCURACY BY AMP 3x SEL TYPES (PSEUDO AMP XVAL & LINEAR FIT)
% KEY: 
% [nj] Model Iterations     ... [scalar, 200]
% [nt] Selectivity Types    ... ["FreqMod","AmpMod","DualMod","InterMod"]
% [nd] Number of Cell Draws ... [various lengths, 5-8 total]
% [nc] Cell Types           ... ["EXC","PV","SOM"]


targetBin = 8; % 8th bin = 50 neuron draws
xVals = 1:5; % pseudo amps
targetSels = [1 3 4]; % freq, dual, inter
figure('Color',[1 1 1],'Theme','Light'); hold on;
for nt = 1 : length(targetSels) % 3 sel types
    for nc = 1 : length(cellType)
        % MDL KEY:  nAmp (5), nIter (200), nSelTypes (3), nDrawBins (5-38), nCellTypes (3)
        yMean = mean(squeeze(PCTall(:,:,targetSels(nt),targetBin,nc))',1).*100;
        ySems = sem(squeeze(PCTall(:,:,targetSels(nt),targetBin,nc))',1).*100;
        xLooped = [xVals flip(xVals)]; 
        yLooped = [yMean+ySems flip(yMean-ySems)]; 
        fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.3,'EdgeColor','none');
        switch nt
            case 1 % freq
                plot(xVals,yMean,'--','LineWidth',1,'Color',colors.(cellType(nc)));  
            case 2 % dual
                plot(xVals,yMean,'-.','LineWidth',1,'Color',colors.(cellType(nc)));  
            case 3 %inter
                plot(xVals,yMean,'-','LineWidth',1,'Color',colors.(cellType(nc)));  
        end
    end
end
% plot([xVals(1) xVals(end)],[20 20],'r--'); %chance performance line
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Amplitude (\mum)','FontName','Arial');
titleTxt = "Frequency decoding performance by amplitude";
title(titleTxt,'FontName','Arial');
subtitle("From all 3 selectivity types, 50 neurons per model")
set(gca,'XLim',[xVals(1) xVals(end)],'YLim',[25 70],'YTick',30:10:80,...
    'XTick',xVals,'XTickLabel',ampValAll);

% FITTING PLOT
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



%% ACCURACY BY AMP 3x SEL TYPES (TRUE AMP XVAL & LOG FIT)
% KEY: 
% [nj] Model Iterations     ... [scalar, 200]
% [nt] Selectivity Types    ... ["FreqMod","AmpMod","DualMod","InterMod"]
% [nd] Number of Cell Draws ... [various lengths, 5-8 total]
% [nc] Cell Types           ... ["EXC","PV","SOM"]


targetBin = 8; % 8th bin = 50 neuron draws
xVals = ampValAll'; % pseudo amps
targetSels = [1 3 4]; % freq, dual, inter
figure('Color',[1 1 1],'Theme','Light'); hold on;

% fake data for correct legend order
plot([-1 -2],[-1 -1],'k-'); plot([1 2],[1 1],'k-.'); % Inter & Dual
plot([-1 -2],[-1 -1],'k--'); % Freq & Amp
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');

for nt = 1 : 3 % only one type here
    for nc = 1 : length(cellType)
        % MDL KEY:  nAmp (5), nIter (200), nSelTypes (3), nDrawBins (5-38), nCellTypes (3)
        yMean = mean(squeeze(PCTall(:,:,targetSels(nt),targetBin,nc))',1).*100;
        ySems = sem(squeeze(PCTall(:,:,targetSels(nt),targetBin,nc))',1).*100;
        xLooped = [xVals flip(xVals)]; 
        yLooped = [yMean+ySems flip(yMean-ySems)]; 
        fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.3,'EdgeColor','none');
        switch nt
            case 1 % freq
                plot(xVals,yMean,'--','LineWidth',1,'Color',colors.(cellType(nc)));  
            case 2 % dual
                plot(xVals,yMean,'-.','LineWidth',1,'Color',colors.(cellType(nc)));  
            case 3 %inter
                plot(xVals,yMean,'-','LineWidth',1,'Color',colors.(cellType(nc)));  
        end
    end
end
% plot([xVals(1) xVals(end)],[20 20],'r--'); %chance performance line
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Amplitude (\mum)','FontName','Arial');
titleTxt = "Frequency decoding performance by amplitude";
title(titleTxt,'FontName','Arial');
subtitle("From all 3 selectivity types, 50 neurons per model")
set(gca,'XLim',[xVals(1) xVals(end)],'YLim',[25 70],'YTick',30:10:80,...
    'XTick',xVals,'XTickLabel',ampValAll);

% % FITTING
% for nc = 1 : length(cellType)
%     yMeans.(cellType(nc)) = mean(squeeze(PCTall2(:,:,nt,maxNeuAvail,nc)),2);
%     [f.(cellType(nc)), gof.(cellType(nc))] = fit(xVals', yMeans.(cellType(nc)),'poly1');
%     plot(xVals',f.(cellType(nc))(xVals).*100,'--','Color',colors.(cellType(nc)));
% end
legend({"Interaction","Dual effect","Freq. effect","PV","SOM","EXC"},'Location','NorthWest','NumColumns',2);



%% PLOT DATA FROM all AMPS (PSEUDO AMP XVAL & LINEAR FIT)
% KEY: 
% [nj] Model Iterations     ... [scalar, 200]
% [nt] Selectivity Types    ... ["FreqMod","AmpMod","DualMod","InterMod"]
% [nd] Number of Cell Draws ... [various lengths, 5-8 total]
% [nc] Cell Types           ... ["EXC","PV","SOM"]


maxDraw = 300;
xVals = 1:5; % pseudo amps
nt = 1; % only one type here
maxNeuAvail = find(squeeze(PCTall2(targetAmp,1,nt,:,3))==0,1)-1; %min cell draw from SOM neurons
figure('Color',[1 1 1],'Theme','Light'); hold on;
for nc = 1 : length(cellType)
    % MDL KEY:  nAmp (5), nIter (200), nCellTypes (3), nCellDraws (5-38), nCellTypes (3)
    yMean = mean(squeeze(PCTall2(:,:,nt,maxNeuAvail,nc))',1).*100;
    ySems = sem(squeeze(PCTall2(:,:,nt,maxNeuAvail,nc))',1).*100;
    xLooped = [xVals flip(xVals)]; 
    yLooped = [yMean+ySems flip(yMean-ySems)]; 
    fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.3,'EdgeColor','none');
    plot(xVals,yMean,'-','LineWidth',1,'Color',colors.(cellType(nc)));    
end
% plot([xVals(1) xVals(end)],[20 20],'r--'); %chance performance line
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Amplitude (\mum)','FontName','Arial');
titleTxt = "Frequency decoding performance by amplitude";
title(titleTxt,'FontName','Arial');
subtitle("From all frequency modulated neurons, 160 per model")
set(gca,'XLim',[xVals(1) xVals(end)],'YLim',[30 80],'YTick',30:10:80,...
    'XTick',xVals,'XTickLabel',ampValAll);

% FITTING
for nc = 1 : length(cellType)
    yMeans.(cellType(nc)) = mean(squeeze(PCTall2(:,:,nt,maxNeuAvail,nc)),2);
    [f.(cellType(nc)), gof.(cellType(nc))] = fit(xVals', yMeans.(cellType(nc)),'poly1');
    plot(xVals',f.(cellType(nc))(xVals).*100,'--','Color',colors.(cellType(nc)));
end
legend({"","EXC","","PV","","SOM","fit","fit","fit"},'Location','NorthWest','NumColumns',2);


%% PLOT DATA FROM all AMPS (TRUE AMP XVAL & LOG FIT)
% KEY: 
% [nj] Model Iterations     ... [scalar, 200]
% [nt] Selectivity Types    ... ["FreqMod","AmpMod","DualMod","InterMod"]
% [nd] Number of Cell Draws ... [various lengths, 5-8 total]
% [nc] Cell Types           ... ["EXC","PV","SOM"]


maxDraw = 300;
xVals = ampValAll'; % pseudo amps
nt = 1; % only one type here
maxNeuAvail = find(squeeze(PCTall2(targetAmp,1,nt,:,3))==0,1)-1; %min cell draw from SOM neurons
figure('Color',[1 1 1],'Theme','Light'); hold on;
for nc = 1 : length(cellType)
    % MDL KEY:  nAmp (5), nIter (200), nCellTypes (3), nCellDraws (5-38), nCellTypes (3)
    yMean = mean(squeeze(PCTall2(:,:,nt,maxNeuAvail,nc))',1).*100;
    ySems = sem(squeeze(PCTall2(:,:,nt,maxNeuAvail,nc))',1).*100;
    xLooped = [xVals flip(xVals)]; 
    yLooped = [yMean+ySems flip(yMean-ySems)]; 
    fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.3,'EdgeColor','none');
    plot(xVals,yMean,'-','LineWidth',1,'Color',colors.(cellType(nc)));    
end
% plot([xVals(1) xVals(end)],[20 20],'r--'); %chance performance line
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Amplitude (\mum)','FontName','Arial');
titleTxt = "Frequency decoding performance by amplitude";
title(titleTxt,'FontName','Arial');
subtitle("From all frequency modulated neurons, 160 per model")
set(gca,'XLim',[xVals(1) xVals(end)],'YLim',[30 80],'YTick',30:10:80,...
    'XTick',xVals,'XTickLabel',ampValAll);

figure('theme','light');
% FITTING
highResX = xVals(1):0.02:xVals(end);
for nc = 1 : length(cellType)
    yMeans.(cellType(nc)) = mean(squeeze(PCTall2(:,:,nt,maxNeuAvail,nc)),2);
    [f.(cellType(nc)), gof.(cellType(nc))] = fit(xVals', yMeans.(cellType(nc)),'log');
    plot(highResX',f.(cellType(nc))(highResX).*100,'--','Color',colors.(cellType(nc)));
end
legend({"","EXC","","PV","","SOM","fit","fit","fit"},'Location','NorthWest','NumColumns',2);


%% PLOT: PERFORMANCE BY SEL TYPE, ALL EXC (1000 ea.)

% KEY: 
% [nj] Model Iterations     ... [scalar, 200]
% [nt] Selectivity Types    ... ["FreqMod","AmpMod","DualMod","InterMod"]
% [nd] Number of Cell Draws ... [various lengths, 5-8 total]
% [nc] Cell Types           ... ["EXC","PV","SOM"]


selType = ["FreqMod","AmpMod","DualMod","InterMod"]; 
newOrder = [4 3 1 2];
cellType = "EXC";
targetAmp = 5;
neuStepSize = 15;
maxDraw = 1000;
figure('Color',[1 1 1],'Theme','Light'); hold on;
for ntn = 1 : length(selType)
    nt = newOrder(ntn);
    for nc = 1 : length(cellType)
        if strcmp(cellType(nc),'EXC')
            maxNeuAvail = maxDraw;
        else % PV | SOM
            maxNeuAvail = sum(is.(cellType(nc)) & is.(selType(nt)));
        end
        nCellDraws = 15:neuStepSize:maxNeuAvail;
        target = cellType(nc);
        nBins = 1:length(nCellDraws);
        yMean = mean(squeeze(PCTall3(targetAmp,:,nt,nBins,nc)),1).*100;
        ySems = sem(squeeze(PCTall3(targetAmp,:,nt,nBins,nc)),1).*100;
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
plot([15 maxDraw],[20 20],'r--');
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Neurons used','FontName','Arial');
title('Frequency decoding performance by selectivity type','FontName','Arial');
subtitle('From random samples of EXC neurons','FontName','Arial');
set(gca,'XLim',[15 1000],'YLim',[10 100],'YTick',10:10:100,'XTick',15:75:1000);
legend({"","Interaction","","Both mod.","","Freq. mod.","","Amp. mod.","chance"},'Location','best');




%% PLOT PC TUNING CURVES FROM MODELS


targetSelType = 4; % "FreqMod","AmpMod","DualMod","InterMod"
targetAmp = 5; % 0.33,0.76,1.6,3.5,7.7 um
visNumPC = 15; % num PCs to visualize
visNumReps = 1; %num reps to visualize
xp = 0.02; % x-padding
numTrials = 12;
% RGB = orderedcolors("gem"); % colorspace for PCs
RGB = mkcolors([0 0 1],[0 1 0],visNumPC);

maxAvailDraws = find(squeeze(PCTall(targetAmp,1,targetSelType,:,targetCellType))==0,1)-1;
targetCellDraw = maxAvailDraws; % cell draw bin index

figure('Theme','light','color',[1 1 1]);
tiledlayout(sqrt(visNumReps),sqrt(visNumReps));
for nr = 1 : visNumReps
    % MDL KEY:  nAmp (5), nIter (200), nCellTypes (3), nCellDraws (5-38), nCellTypes (3)
    mdl = MDLall{targetAmp,nr,targetSelType,targetCellDraw,targetCellType}; % nIter (200), nCellTypes (4), nCellDraws (5-38), nCellTypes (3)
    PCresps = zeros(numTrials,nFreq,visNumPC);
    nexttile(nr); hold on;
    for np = 1 : visNumPC
        xpl = xp+(np-(round(visNumPC/2)+1))*xp; % local iteration
        for nf = 1 : nFreq
            PCresps(:,nf,np) = mdl.X(mdl.Y == freqVal(nf),np);
            ySem = sem(PCresps(:,nf,np));
            yMean = mean(PCresps(:,nf,np));
            plot([nf nf]+xpl,[yMean-ySem yMean+ySem],'-','Color',RGB(np,:)); %sem lines
        end
        plot((1:5)+xpl, mean(squeeze(PCresps(:,:,np)),1),'-','Color',RGB(np,:));%mean line
    end
    set(gca,'XLim',[0.5 5.5],'XTick',1:5,'XTickLabel',freqVal);
    ylabel('z-scored "response"','FontName','Arial');
    xlabel('frequency (Hz)','FontName','Arial');
    title(strcat(cellType(targetCellType)," model rep #",num2str(nr)),'FontName','Arial');
    subtitle("PC1 (Blue) ... PC15 (Green)");
    hold off;
end


%% hOW DOES PC TUNING STRENGTH CHANGE AS MORE NEURONS ARE ADDED

targetSelTypes = [1 3 4]; % freqSel, dualSel, interaction
tAmp = 5;
tRep = 1;
numTrials = 12; % minTrialNum * trainPortion
visNumPC = 3;
xVals = 15:5:120;
figure('Theme','light','color',[1 1 1]); 
tiledlayout(3,3);

for nt = 1 : length(targetSelTypes) % LOOP: selectivity type 
	for nc = 1 : size(MDLall,5)  % LOOP: cell type
        nexttile; hold on;
        for np = 1 : visNumPC % LOOP: num PCs
            selVals = zeros(size(MDLall,2),size(MDLall,4)); %200 x 22
		    for nb = 1 : size(MDLall,4) % LOOP: num draw bins
			    for nr = 1 : size(MDLall,2) % LOOP: num reps
				    mdl = MDLall{tAmp,nr,targetSelTypes(nt),nb,nc};
                    if isempty(mdl)
                        break
                    end
                    PCresps = zeros(numTrials,nFreq);
        			for nf = 1 : nFreq 
            			PCresps(:,nf) = mdl.X(mdl.Y == freqVal(nf),np);
                    end
                    selVals(nr,nb) = quickSelectivity(squeeze(mean(PCresps,1)));
                end
                if isempty(mdl)
                    break
                end
            end
            usableCols = ~all(selVals == 0,1);
            meanSelVals = mean(selVals(:,usableCols),1);
            switch np
                case 1
                    plot(xVals(usableCols),meanSelVals,'-','Color',colors.(cellType(nc)));
                case 2
                    plot(xVals(usableCols),meanSelVals,'--','Color',colors.(cellType(nc)));
                case 3
                    plot(xVals(usableCols),meanSelVals,':','Color',colors.(cellType(nc)));
            end
            set(gca,'XLim',[min(xVals(usableCols))-2 max(xVals(usableCols))+2],...
                'XTick',[15:15:max(xVals(usableCols))])
        end
        if nt == 2 && nc == 1
            ylabel("PC Selectivity Strength","FontName","Arial");
        end
        if nt == 3 && nc == 2
            xlabel("Number of neurons used","FontName","Arial");
        end
        if nc == 2
            title(strcat("Neurons with ",selType(targetSelTypes(nt))),"FontName","Arial");
        end
    end
end


%% hOW DOES PC VARIANCE PER FREQ CHANGE AS MORE NEURONS ARE ADDED

targetSelTypes = [1 3 4]; % freqSel, dualSel, interaction
tAmp = 5;
tRep = 1;
numTrials = 12; % minTrialNum * trainPortion
visNumPC = 3;
xVals = 15:5:120;
figure('Theme','light','color',[1 1 1]); 
tiledlayout(3,3);

for nt = 1 : length(targetSelTypes) % LOOP: selectivity type 
	for nc = 1 : size(MDLall,5)  % LOOP: cell type
        nexttile; hold on;
        for np = 1 : visNumPC % LOOP: num PCs
            varVals = zeros(size(MDLall,2),size(MDLall,4)); %200 x 22
		    for nb = 1 : size(MDLall,4) % LOOP: num draw bins
			    for nr = 1 : size(MDLall,2) % LOOP: num reps
				    mdl = MDLall{tAmp,nr,targetSelTypes(nt),nb,nc};
                    if isempty(mdl)
                        break
                    end
                    PCresps = zeros(numTrials,nFreq);
        			for nf = 1 : nFreq 
            			PCresps(:,nf) = mdl.X(mdl.Y == freqVal(nf),np);
                    end
                    varVals(nr,nb) = mean(std(PCresps,1));
                end
                if isempty(mdl)
                    break
                end
            end
            usableCols = ~all(varVals == 0,1);
            meanVarVals = mean(varVals(:,usableCols),1);
            switch np
                case 1
                    plot(xVals(usableCols),meanVarVals,'-','Color',colors.(cellType(nc)));
                case 2
                    plot(xVals(usableCols),meanVarVals,'--','Color',colors.(cellType(nc)));
                case 3
                    plot(xVals(usableCols),meanVarVals,':','Color',colors.(cellType(nc)));
            end
            set(gca,'XLim',[min(xVals(usableCols))-2 max(xVals(usableCols))+2],...
                'XTick',[15:15:max(xVals(usableCols))])
        end
        if nt == 2 && nc == 1
            ylabel("PC tuning variation (avg. std)","FontName","Arial");
        end
        if nt == 3 && nc == 2
            xlabel("Number of neurons used","FontName","Arial");
        end
        if nc == 2
            title(strcat("Neurons with ",selType(targetSelTypes(nt))),"FontName","Arial");
        end
    end
end


%% EVALUATE PC CONTRIBUTIONS TO MODEL

% SELECT EXAMPLE MODEL
targetAmp = 5; % 0.33,0.76,1.6,3.5,7.7 um
targetSelType = 4; % "FreqMod","AmpMod","DualMod","InterMod"
targetCellType = 2; % [EXC, PV, SOM]
targetRep = 3; % (1-200)
% Find model with largest number of neurons included
% maxAvailDraws = find(squeeze(PCTall(targetAmp,1,targetSelType,:,targetCellType))==0,1)-1;
% targetCellDraw = maxAvailDraws; % cell draw bin index
targetCellDraw = 8;
mdl = MDLall{targetAmp,targetRep,targetSelType,targetCellDraw,targetCellType};

% Verify parameter sizes
numClasses = numel(unique(mdl.Y)); % number of classes
numPreds = size(mdl.Mu, 2); % number of predictors

% Extract pairwise discriminant coefficients 
coefMat = zeros(numClasses, numClasses, numPreds);
for i = 1:numClasses
    for j = 1:numClasses
        if ~isempty(mdl.Coeffs(i,j).Linear)
            coefMat(i,j,:) = mdl.Coeffs(i,j).Linear;
        end
    end
end

% Average absolute coefficient magnitude across all pairwise discriminants
coefImportance = squeeze(mean(abs(coefMat), [1 2], 'omitnan'));

% Normalize both metrics for comparability
coefImportance = coefImportance / max(coefImportance);
deltaImportance = mdl.DeltaPredictor(:) / max(mdl.DeltaPredictor);

% ---- 2. Plot them together ----
figure;
hold on
b1 = bar(1:numPreds, coefImportance, 0.45, 'FaceColor', [0.3 0.6 0.9], 'DisplayName', '|Coefs| (boundary strength)');
b2 = bar((1:numPreds)+0.45, deltaImportance, 0.45, 'FaceColor', [0.9 0.5 0.2], 'DisplayName', 'ΔPredictor (performance drop)');
hold off
xlabel('Predictor (Neuron)','FontName','Arial');
ylabel('Normalized Importance','FontName','Arial');
title('Predictor Impact: Coefficient vs ΔPredictor Measures','FontName','Arial');
legend('Location', 'best');



%% BIAS OVER/UNDER PREDICTION ANALYSIS 

pValThresh = 0.05;

% cd('D:\Tactile_Synchrony\2P-Data\Cleaned\datasets');
load('modelResults_dataset2_ampDecode_predictionAnalysis.mat') % 'PCTall','predFreqAll','testFreqAll'
% dimension IDs: nFreq, nIter, selType, 1, cellType

% FOR Prediction analysis with all selective types pooled, re-run first
% section code above with selType as 'allMod'

targSel = 1;
xSp = [-0.2 0 0.2]; % x-axis plot spacing
allDiffs = zeros(size(PCTall,1),size(PCTall,2),size(PCTall,3),size(PCTall,5)); % 5 x 200 x 4 x 3

figure('Color',[1 1 1]); hold on;
plot([0 6],[0 0],'k:','LineWidth',0.2);
for na = 1 : size(PCTall,1) % loop amplitude (5)
    for nc = 1 : size(PCTall,5) % loop Cell Type (3)
        for ni = 1 : size(PCTall,2) % loop nIter (200)
            locPredVals = predFreqAll{na,ni,targSel,1,nc}; % predicted amplitudes
            locTrueVals = testFreqAll{na,ni,targSel,1,nc}; % true amplitudes
            locDiffs = locPredVals - locTrueVals; % differences (predicted - true)
            locDiffs(locDiffs == 0) = []; % remove correct values (0)
            allDiffs(na,ni,targSel,nc) = mean(locDiffs);
        end
        locX = na + xSp(nc); % x-position
        locY = mean(squeeze(allDiffs(na,:,targSel,nc)),'omitnan'); % average difference
        locCI = mkCI(squeeze(allDiffs(na,:,targSel,nc))); % 95 percent condifidence intervals
        errorbar(locX,locY,locY-locCI(1),'vertical','Color',colors.(cellType(nc)),'LineWidth',1);
        plot(locX,locY,'.','Color',colors.(cellType(nc)),'MarkerSize',10);
    end
end
set(gca,'XTick',1:5,'XTickLabels',ampValAll,'XLim',[0.6 5.4]);
xlabel('Amplitude (\mum)','FontName','Arial');
ylabel({'Average prediction error (Hz)','(predicted - true)'},'FontName','Arial');
title('Frequency decoding bias','FontName','Arial');

% EVALUATE WITH MODEL FITS AND PLOT
for nc = 1 : length(cellType)
    disp(cellType(nc));
    locYs = squeeze(allDiffs(:,:,targSel,nc)); % average difference values, amp by reps
    locY_vect = locYs(~isnan(locYs)); % remove NaN and vectorize
    locXs = repmat((1:5)',1,size(allDiffs,2)); % average difference values, amp by reps
    locX_vect = locXs(~isnan(locYs)); % remove NaN and vectorize

    % LINEAR FIT
    locMdl = fitlm(locX_vect,locY_vect); % linear model
    AIC_linear = locMdl.ModelCriterion.AIC; % AIC
    fitSigLinear = locMdl.ModelFitVsNullModel.Pvalue < pValThresh; % significance of linear fit
    slope_Linear = locMdl.Coefficients.Estimate(2);
    rsq_Linear = locMdl.Rsquared.Adjusted;
    pval_Linear = locMdl.ModelFitVsNullModel.Pvalue;

    % POWER 2 FIT
    [fitobject,gof,output] = fit(locX_vect,locY_vect,'power2'); % power2 model
    AIC_power = output.numobs * log(gof.sse / output.numobs) + 2 * output.numparam; % AIC
    confBounds = confint(fitobject); % 95% confidence bounds of coefficients
    fitSigPower = all(confBounds(1,:).*confBounds(2,:) > 0); % check for same sign in coeff bounds

    % COMPARE MODELS
    xFine = linspace(min(locX_vect),max(locX_vect),100)';
    if fitSigLinear == true && fitSigPower == false
        % ONLY Linear model significant
        yFit = predict(locMdl,xFine);
        plot(xFine+xSp(nc),yFit,'--','Color',colors.(cellType(nc)));
    elseif fitSigLinear == false && fitSigPower == true 
        % ONLY Power model significant
        yFit = feval(fitobject, xFine); 
        yPI = predint(fitobject, xFine, 0.95); % 95% prediction intervals
        coefBounds = confint(fitobject, 0.95); % coeff. conf. bounds
        plot(xFine+xSp(nc),yFit,'--','Color',colors.(cellType(nc)));
    elseif fitSigLinear == true && fitSigPower == true 
        % BOTH models are significant, compare AIC
        if AIC_linear > AIC_power 
            % Linear model has higher AIC
            yFit = predict(locMdl,xFine);
            plot(xFine+xSp(nc),yFit,'--','Color',colors.(cellType(nc)));
        elseif AIC_linear < AIC_power 
            % Power model has higher AIC
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

