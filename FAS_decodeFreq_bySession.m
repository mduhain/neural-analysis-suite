% FAS_decodeFreq_bySession.m
%
% within-session frequency decoding script
% outputs saved and then used in: FAS_neuCorrByFreqDecoding.m


% USER INPUTS (CHANGE THESE)
nIter = 200; % Number of model iterations to run
trainPortion = 0.8; % Percentage of trials to be used in training (remaining for testing)
usePCA = true; % Use PCA step prior to model training
kStart = 1; % First PC to include
kComp = 15; % Total number of PCs to use in model training / testing
neuStepSize = 15; % Num neurons to add to model each step
gammaSteps = 0; % Gamma value(s) for LDA
deltaSteps = 0; % Delta values(s) for LDA; [0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1];
targetAmp = 5; % Amplitude value to decode with; [0.33 0.76 1.6 3.5 7.7] um
onlyEXC = false; % Decoding with only EXC neurons in sessions
lateRespWin = false; % Decoding using the later response window (3-4 sec post-stim)
onlyFreqSel = false; % Decoding with only frequency-selective units

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% CREATE OUTPUT ARRAYS
PCTall = zeros(nIter,length(sessIDs),length(gammaSteps),length(deltaSteps),2); % (nj,ns,ng,ne,sess vs pseudopop)
MDLall = cell(nIter,length(sessIDs),length(gammaSteps),length(deltaSteps),2); % (nj,ns,ng,nd)
CVEall = zeros(nIter,length(sessIDs),length(gammaSteps),length(deltaSteps),2);
PCFall = cell(nIter,length(sessIDs),length(gammaSteps),length(deltaSteps),2); % (nj,ns,ng,nd)
PredAll = cell(nIter,length(sessIDs),length(gammaSteps),length(deltaSteps),2); % all predicted frequencies
TrueAll = cell(nIter,length(sessIDs),length(gammaSteps),length(deltaSteps),2); % all true frequencies

for ns = 1:length(sessIDs) %nt = 1 : length(selType)
    locSess = char(sessIDs(ns)); % local session name 'm000_00'
    locIdx = is.mouseNum(:,strcmp(mouseNums,locSess(1:4))); % all neurons/rows from this mouse
    disp(strcat("Session ",num2str(ns),": ",locSess));
    % disp(selType(nt));

    % nCellDraws = kComp : neuStepSize : maxNeuAvail;
    trainTrial = round(minTrialNum*trainPortion);
    testTrial = round(minTrialNum*(1-trainPortion));
    trainFreq = repmat(freqVal,trainTrial,1);
    testFreq = repmat(freqVal,testTrial,1);
    
    % loc cell types within session
    % ! PV or SOM, one of these arrays is always empty, sometimes EXC is empty

    % DATASET 1 FREQUENCY SELECTIVE NEURONS
    if multiAmp == true && onlyFreqSel == true
        locPVs = find(is.sessID(:,ns) & is.PV & is.allMod); % loc neurons from this session
        locSOMs = find(is.sessID(:,ns) & is.SOM & is.allMod);
        locEXCs = find(is.sessID(:,ns) & is.EXC & is.allMod);
        locNumINH = length(locPVs) + length(locSOMs);
        otherPVs = find(locIdx & is.PV & is.allMod); % loc neurons from this mouse (but not this session)
        otherSOMs = find(locIdx & is.SOM & is.allMod);
        otherEXCs = find(locIdx & is.EXC  & is.allMod);
    % DATASET 2 FREQUENCY SELECTIVE NEURONS
    elseif multiAmp == false && onlyFreqSel == true
        locPVs = find(is.sessID(:,ns) & is.PV & is.selective); % loc neurons from this session
        locSOMs = find(is.sessID(:,ns) & is.SOM & is.selective);
        locEXCs = find(is.sessID(:,ns) & is.EXC & is.selective);
        locNumINH = length(locPVs) + length(locSOMs);
        otherPVs = find(locIdx & is.PV & is.selective); % loc neurons from this mouse (but not this session)
        otherSOMs = find(locIdx & is.SOM & is.selective);
        otherEXCs = find(locIdx & is.EXC & is.selective);
    % DATASET 1 ALL NEURONS
    elseif multiAmp == true && onlyFreqSel == false
        locPVs = find(is.sessID(:,ns) & is.PV); % loc neurons from this session
        locSOMs = find(is.sessID(:,ns) & is.SOM);
        locEXCs = find(is.sessID(:,ns) & is.EXC);
        locNumINH = length(locPVs) + length(locSOMs);
        otherPVs = find(locIdx & is.PV); % loc neurons from this mouse (but not this session)
        otherSOMs = find(locIdx & is.SOM); 
        otherEXCs = find(locIdx & is.EXC); 
    % DATASET 2 ALL NEURONS
    elseif multiAmp == false && onlyFreqSel == false 
        locPVs = find(is.sessID(:,ns) & is.PV); % loc neurons from this session
        locSOMs = find(is.sessID(:,ns) & is.SOM); 
        locEXCs = find(is.sessID(:,ns) & is.EXC); 
        locNumINH = length(locPVs) + length(locSOMs);
        otherPVs = find(locIdx & is.PV ); % loc neurons from this mouse (but not this session)
        otherSOMs = find(locIdx & is.SOM);
        otherEXCs = find(locIdx & is.EXC);
    end

    % count of neurons within this session
    if length(locEXCs) >= locNumINH
        numNeuSess = length(locEXCs);
    elseif isempty(locEXCs) % cre-dep. gcamp expression INH only
        numNeuSess = locNumINH;
        if onlyEXC == true
            disp('No EXC neurons this sess.')
            continue
        end
    elseif length(locEXCs) <= locNumINH
        numNeuSess = locNumINH;
    else
        error('un-expected condition met');
    end

    if usePCA == true && numNeuSess < kComp
        disp('fewer freq.sel. neurons than min PC number');
        continue;
    end
    
    % reiterate random draws
    for nj = 1:nIter
        
        % build array: session neurons (all EXC) or (some EXC + all INH)
        if onlyEXC == true
            sessNeus = locEXCs; % use all EXC neurons available
        elseif onlyEXC == false
            locExcNum = length(locEXCs) - locNumINH; % use fewer EXC neurons to make room for INH neurons
            subSampledEXC = randsample(locEXCs,locExcNum,false); 
            sessNeus = cat(1,locPVs,locSOMs,subSampledEXC); % combine INH and EXC neurons 
        end

        % build array: pseudo populations to match sess cells
        selectedPVs = randsample(otherPVs, length(locPVs), false); % random draw from total population
        selectedSOMs = randsample(otherSOMs, length(locSOMs), false); % random draw from total population
        selectedEXCs = randsample(otherEXCs, length(locEXCs), false); % random draw from total population
        if onlyEXC == true
            selectedCells = selectedEXCs; % ONLY USE EXC NEURONS FOR DECODING
        elseif onlyEXC == false
            selectedCells = cat(1,selectedPVs,selectedSOMs,selectedEXCs); % random sampled pseudo-pop
        end
       
        % pre-allocate training and testing data array
        trainRespS = zeros(nFreq,trainTrial,numNeuSess); % session
        testRespS = zeros(nFreq,testTrial,numNeuSess);
        trainRespP = zeros(nFreq,trainTrial,numNeuSess); % pseudo pop
        testRespP = zeros(nFreq,testTrial,numNeuSess);

        % loop over selected cells
        for nk = 1:length(sessNeus)

            % FROM THIS SESSION
            if lateRespWin == true
                cellResp = Tc.lateResps{sessNeus(nk)}; % 3-4 SEC WINDOW
            elseif lateRespWin == false
                cellResp = Tc.responses{sessNeus(nk)}; % 0-1 SEC WINDOW
            end
            
            cellFreq = Tc.frequencies{sessNeus(nk)};
            if multiAmp
                cellAmp = Tc.amplitudes{sessNeus(nk)};
            end
                      
            for nl = 1:nFreq
                % randomly draw trials at each frequency, no replacement
                if multiAmp
                    trialNum = find(cellFreq == freqVal(nl) & cellAmp == ampValAll(targetAmp));
                else
                    trialNum = find(cellFreq == freqVal(nl));
                end
                % populate training and testing arrays
                selectedTrials = randsample(trialNum, trainTrial+testTrial, false);
                trainRespS(nl,:,nk) = cellResp(selectedTrials(1:trainTrial));
                testRespS(nl,:,nk) = cellResp(selectedTrials(trainTrial+1:end));
            end

            % FROM PSEUDO-POP
            if lateRespWin == true
                cellResp = Tc.lateResps{selectedCells(nk)}; % 3-4 SEC WINDOW
            elseif lateRespWin == false
                cellResp = Tc.responses{selectedCells(nk)}; % 0-1 SEC WINDOW
            end
            

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
                % populate training and testing arrays
                selectedTrials = randsample(trialNum, trainTrial+testTrial, false);
                trainRespP(nl,:,nk) = cellResp(selectedTrials(1:trainTrial));
                testRespP(nl,:,nk) = cellResp(selectedTrials(trainTrial+1:end));
            end
        end

        % Reshape matricies
        trainRespS = reshape(trainRespS, nFreq*trainTrial, numNeuSess);
        testRespS = reshape(testRespS, nFreq*testTrial, numNeuSess);
        trainRespP = reshape(trainRespP, nFreq*trainTrial, numNeuSess);
        testRespP = reshape(testRespP, nFreq*testTrial, numNeuSess);

        % Z-Scoring 
        for nz = 1 : numNeuSess
            % session
            mu = mean(trainRespS(:,nz));
            sigma = std(trainRespS(:,nz));
            trainRespS(:,nz) = (trainRespS(:,nz) - mu) ./ sigma;
            testRespS(:,nz) = (testRespS(:,nz) - mu) ./ sigma;
            % pseudo-pop
            mu = mean(trainRespP(:,nz));
            sigma = std(trainRespP(:,nz));
            trainRespP(:,nz) = (trainRespP(:,nz) - mu) ./ sigma;
            testRespP(:,nz) = (testRespP(:,nz) - mu) ./ sigma;
        end

        % Percent variance explained
        [coeffS,scoreS,latentS] = pca(trainRespS);
        % PC.(cellType(nc)).coeff{ni} = coeff;
        % PC.(cellType(nc)).score{ni} = score;
        % PC.(cellType(nc)).latent{ni} = latent;
        explainedS = 100 * latentS / sum(latentS);
        cumFracS = cumsum(sort(explainedS,'descend'))/sum(explainedS);
        % CVEall(nj,ns,ng,ne) = sum(explained(1:15));
        CVEall(nj,ns,1,1,1) = find(cumFracS >= 0.80, 1); %num of PCs to explain 80% of variance (sess)
        % --------------------
        [coeffP,scoreP,latentP] = pca(trainRespP);
        explainedP = 100 * latentP / sum(latentP);
        cumFracP = cumsum(sort(explainedP,'descend'))/sum(explainedP);
        CVEall(nj,ns,1,1,2) = find(cumFracP >= 0.80, 1); %num of PCs to explain 80% of variance (pseudopop)
        % ---------------------
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

        % TRAIN WITH PCA
        if usePCA == true
            % Session data
            trainRespS_red = scoreS(:,kStart:kComp);
            PCtunings = zeros(length(freqVal),kComp);
            PCvariation = zeros(length(freqVal),kComp);
            for nf = 1 : length(freqVal)
                PCtunings(nf,kStart:kComp) = mean(trainRespS_red(trainFreq == freqVal(nf),:),1);
                PCvariation(nf,kStart:kComp) = std(trainRespS_red(trainFreq == freqVal(nf),:),1);
            end
            for ng = 1 : length(gammaSteps)
                for ne = 1 : length(deltaSteps)
                    mdl = fitcdiscr(trainRespS_red,trainFreq,'DiscrimType','linear','Gamma',gammaSteps(ng),...
                        'Delta',deltaSteps(ne));
                    predFreq = predict(mdl,testRespS*coeffS(:,kStart:kComp));
                    % collect output
                    PredAll{nj,ns,ng,ne,1} = predFreq;
                    TrueAll{nj,ns,ng,ne,1} = testFreq;
                    PCTall(nj,ns,ng,ne,1) = mean(predFreq == testFreq); % accuracy
                    MDLall{nj,ns,ng,ne,1} = mdl; % model parameters
                    % PCFall{nj,ns,ng,ne,1} = cat(3, PCtunings, PCvariation); % accuracy
                    PCFall{nj,ns,ng,ne,1} = coeffS; % coeffs from PCA
                end
            end
            % PSeudo-pop data
            trainRespP_red = scoreP(:,kStart:kComp);
            PCtunings = zeros(length(freqVal),kComp);
            PCvariation = zeros(length(freqVal),kComp);
            for nf = 1 : length(freqVal)
                PCtunings(nf,kStart:kComp) = mean(trainRespP_red(trainFreq == freqVal(nf),:),1);
                PCvariation(nf,kStart:kComp) = std(trainRespP_red(trainFreq == freqVal(nf),:),1);
            end
            for ng = 1 : length(gammaSteps)
                for ne = 1 : length(deltaSteps)
                    mdl = fitcdiscr(trainRespP_red,trainFreq,'DiscrimType','linear','Gamma',gammaSteps(ng),...
                        'Delta',deltaSteps(ne));
                    predFreq = predict(mdl,testRespP*coeffP(:,kStart:kComp));
                    % collect output
                    PCTall(nj,ns,ng,ne,2) = mean(predFreq == testFreq); % accuracy
                    MDLall{nj,ns,ng,ne,2} = mdl; % model parameters
                    % PCFall{nj,ns,ng,ne,2} = cat(3, PCtunings, PCvariation); % accuracy
                    PCFall{nj,ns,ng,ne,2} = coeffS;
                end
            end

        % TRAIN WITHOUT PCA
        elseif usePCA == false
            for ng = 1 : length(gammaSteps)
                for ne = 1 : length(deltaSteps)
                    % Session
                    mdl = fitcdiscr(trainRespS,trainFreq,'DiscrimType','linear','Gamma',gammaSteps(ng),...
                        'Delta',deltaSteps(ne));
                    predFreq = predict(mdl,testRespS);
                    PredAll{nj,ns,ng,ne,1} = predFreq;
                    TrueAll{nj,ns,ng,ne,1} = testFreq;
                    PCTall(nj,ns,ng,ne,1) = mean(predFreq == testFreq); % accuracy
                    MDLall{nj,ns,ng,ne,1} = mdl; % model parameters
                    % Pseudo-Pop
                    mdl = fitcdiscr(trainRespP,trainFreq,'DiscrimType','linear','Gamma',gammaSteps(ng),...
                        'Delta',deltaSteps(ne));
                    predFreq = predict(mdl,testRespP);
                    PCTall(nj,ns,ng,ne,2) = mean(predFreq == testFreq); % accuracy
                    MDLall{nj,ns,ng,ne,2} = mdl; % model parameters
                end
            end
        end
    end 
    % end
    fprintf(newline)
end


allPVals = zeros(ns,1);
allOS = zeros(ns,1); % observed stats
for n = 1 : size(PCTall,2)
    array1 = PCTall(:,n,1,1,1);
    array2 = PCTall(:,n,1,1,2);
    [pVal, obsStat, permStats] = permTest2sample(array1, array2, 5000);
    allPVals(n) = pVal;
    allOS(n) = obsStat;
end

PCT = squeeze(mean(PCTall,1));
sig1 = allPVals < pValThresh & allOS < 0;
sig2 = allPVals < pValThresh & allOS > 0;

% if usePCA
%     wPCA = struct();
%     wPCA.PCTall = PCTall;
%     wPCA.MDLall = MDLall;
%     wPCA.allPVals = allPVals;
%     wPCA.allOS = allOS;
% else
%     woPCA = struct();
%     woPCA.PCTall = PCTall;
%     woPCA.MDLall = MDLall;
%     woPCA.allPVals = allPVals;
%     woPCA.allOS = allOS;
% end

% save out
if onlyEXC == true && lateRespWin == true
    PCTallEXC_late = squeeze(PCTall);
elseif onlyEXC == false && lateRespWin == true   
    PCTall_late = squeeze(PCTall);
elseif onlyEXC == true && lateRespWin == false
    PCTallEXC_baseline = squeeze(PCTall);
elseif onlyEXC == false && lateRespWin == false 
    PCTall_baseline = squeeze(PCTall);
end

%% PV+EXC vs. SOM+EXC
figure('Color',[1 1 1]); hold on;
% tiledlayout(4,4);

% nexttile([3 4]); hold on;
accPV = [];
accSOM = [];
for ns = 1 : length(sessIDs)
    if PCT(ns) ~= 0 % session with non-zero avg accuracy
        if any(is.sessID(:,ns) & is.PV) % PV + EXC session
            accPV = cat(1,accPV,PCT(ns,1)*100);
        elseif any(is.sessID(:,ns) & is.SOM) % SOM + EXC session
            accSOM = cat(1,accSOM,PCT(ns,1)*100);
        end
    end
end
histogram(accPV,'Normalization','percentage','NumBins',10,'FaceAlpha',0.5,'FaceColor',colors.PV);
histogram(accSOM,'Normalization','percentage','NumBins',10,'FaceAlpha',0.5,'FaceColor',colors.SOM);
ylabel('% of sessions','FontName','Arial');
xlabel('Decoding accuracy (%)','FontName','Arial');
title('Within-session frequency decoding','FontName','Arial');
set(gca,'XLim',[7 77],'YLim',[0 21]);

% nexttile([1 4]); hold on;
meanPV = mean(accPV);
ciPV = mkCI(accPV);
errorbar(meanPV,20,ciPV(2)-meanPV,'horizontal','Color',colors.PV,'LineWidth',1);
plot(meanPV,20,'.','MarkerSize',10,'Color',colors.PV);
meanSOM = mean(accSOM);
ciSOM = mkCI(accSOM);
errorbar(meanSOM,20,ciSOM(2)-meanSOM,'horizontal','Color',colors.SOM,'LineWidth',1);
plot(meanSOM,20,'.','MarkerSize',10,'Color',colors.SOM);

legend({'PV+EXC','SOM+EXC',''},'Box','off','FontName','Arial');


%% Confusion Matrix (old version)

% Load in Data
load('PCT_dataset1_withinSess_freqSel.mat','PredAll','TrueAll');
PredFS = PredAll;
TrueFS = TrueAll;
load('PCT_dataset1_withinSess_allAvail.mat','PredAll','TrueAll');
PredAA = PredAll;
TrueAA = TrueAll;

% Plotting
figure;
% tiledlayout(1,2);

% % FREQUENCY-SELECTIVE UNITS
% nexttile();
% allPredValsFS = [];
% allTrueValsFS = [];
% for ns = 1:length(sessIDs)
%     targSess = ns;
%     % check if session was fully trained / predicted with
%     if isempty(PredFS{1,targSess,1,1,1})
%         continue
%     end
%     allPredValsFS = cat(1,allPredValsFS,PredFS{:,targSess,1,1,1});
%     allTrueValsFS = cat(1,allTrueValsFS,TrueFS{:,targSess,1,1,1});
% end
% confusionchart(allTrueValsFS,allPredValsFS,'Normalization','total-normalized',...
%     "DiagonalColor",'b','OffDiagonalColor','r');
% xlabel('Predicted Frequency (Hz)');
% ylabel('True Frequency (Hz)');

% ALL AVAILABLE UNITS
% nexttile();
allPredValsAA = [];
allTrueValsAA = [];
for ns = 1:length(sessIDs)
    targSess = ns;
    % check if session was fully trained / predicted with
    % using FS data, since there are more AA neurons, more fully trained
    % sessions which dont have a paired FS model set.
    if isempty(PredFS{1,targSess,1,1,1})
        continue
    end
    allPredValsAA = cat(1,allPredValsAA,PredAA{:,targSess,1,1,1});
    allTrueValsAA = cat(1,allTrueValsAA,TrueAA{:,targSess,1,1,1});
end
confusionchart(allTrueValsAA,allPredValsAA,'Normalization','column-normalized',...
    "DiagonalColor",'k','OffDiagonalColor','k');
xlabel('Predicted Frequency (Hz)');
ylabel('True Frequency (Hz)');


%% new alternative confusion chart (2026-07-21)

% compute confusion matrix and class order
[C,order] = confusionmat(allTrueValsAA,allPredValsAA);

% column-normalized matrix (columns sum to 1)
colSums = sum(C,1);
colNorm = C ./ colSums;

% total-normalized matrix (each entry divided by total sum)
totalNorm = C / sum(C(:));

% display
figure('Color',[1 1 1]);
imagesc(colNorm.*100);

% % Custom grayscale colormap from [0.9 0.9 0.9] to [0.1 0.1 0.1] with CLim [0 70]
% nC = 256; % number of colors
% cmap = [linspace(0.9,0.1,nC)' repmat(linspace(0.9,0.1,nC)',1,2)];

colormap(parula);
clim([0 70]); % set color limits
colorbar; % optional: show colorbar

% adjust
set(gca,'YTick',1:5,'YTickLabel',order,'XTickLabel',order,'FontName','Arial');
xlabel('Predicted Frequency (Hz)','FontName','Arial');
ylabel('True Frequency (Hz)','FontName','Arial');
title('Average within-session frequency decoding','FontName','Arial');
subtitle('Column-normalized (color), total-normalized (label)','FontName','Arial');
for nv = 1 : 5
    locVal = strcat(num2str(totalNorm(nv,nv)*100,3),"%");
    text(nv,nv,locVal,'HorizontalAlignment','center','VerticalAlignment','middle')
end



%% Dataset 1: Within-Session Decoding, frequency-selective vs. all-available

% load('PCT_dataset1_withinSess_freqSelVsAllAvail.mat'); % 'PCT_FS' & 'PCT_AA'
% PCT_FS = PCT_FS(:,1) * 100; % remove pseudo-pop data
% PCT_AA = PCT_AA(:,1) * 100; % remove pseudo-pop data


% 2026-07-14 (validation run)
load('PCT_dataset1_withinSess_freqSel.mat','PCTall');
PCT_FS = mean(PCTall(:,:,1,1,1),1)' * 100;
load('PCT_dataset1_withinSess_allAvail.mat','PCTall');
PCT_AA = mean(PCTall(:,:,1,1,1),1)' * 100;


idx = (PCT_FS(:) ~= 0) & (PCT_AA(:) ~= 0);

% TWO BAR PLOTS, FS vs AA
figure('Color',[1 1 1]); hold on;

bar(1,mean(PCT_FS(idx)),'FaceColor',[0.75 0.75 0.75]);
locCIfs = mkCI(PCT_FS(idx));
errorbar(1,mean(PCT_FS(idx)),mean(PCT_FS(idx))-locCIfs(1),'Color',[0 0 0]);
%swarmchart(ones(length(PCT_FS(idx)),1),PCT_FS(idx),'filled','o','MarkerFaceColor',[0 0 0],'SizeData',4);

bar(2,mean(PCT_AA(idx)),'FaceColor',[0.9 0.9 0.9]);
locCIaa = mkCI(PCT_AA(idx));
errorbar(2,mean(PCT_AA(idx)),mean(PCT_AA(idx))-locCIaa(1),'Color',[0 0 0]);
%swarmchart(ones(length(PCT_AA(idx)),1)+1,PCT_AA(idx),'filled','o','MarkerFaceColor',[0 0 0],'SizeData',4);

set(gca,'XTick',[],'XLim',[0.2 2.8],'YLim',[35 60]);
ylabel('Frequency decoding accuracy (%)','FontName','Arial');
title("Model performance within sessions",'FontName','Arial');

% 2-sided Kolmogorov-Smirnov test
[h,p,ks2stat] = kstest2(PCT_FS(idx),PCT_AA(idx));
yv = max([locCIfs; locCIaa]) + 2; % y-value for figure annotation
if p < 0.05
    plot([1 2],[yv yv],'k-');
    text(1.5,yv,'*','FontName','Arial','HorizontalAlignment','center','VerticalAlignment','bottom');
else
    plot([1 2],[yv yv],'k-');
    text(1.5,yv,'n.s.','FontName','Arial','HorizontalAlignment','center','VerticalAlignment','bottom');
end
legend({'Frequency-Selective','','All Available',''},'Box','off');


%% PLOTTING

% SCATTER (SESS vs. PSEUDO) + 2 EXAMPLES
% Do indv. sessions outpreform pseudo-populations in freq. decoding?
figure('Color',[1 1 1]);
% plot(PCT(:,1)*100,allOS*100,'ko');
% ylabel('change in decoding accuracy % (sess v. pseudo)','FontName','Arial');
subplot(2,3,[1 2 4 5]); hold on;
plot([0 100],[0 100],'k:');

% PCT ratio (sess:pseudo)
PCTr = PCT(:,1) ./ PCT(:,2);
histogram(PCTr(~sig1 & ~sig2),'BinEdges',0:0.02:1,'Normalization','percentage',...
    'FaceAlpha',0.5,'FaceColor',[.2,.2,.2],'EdgeAlpha',0.5,'EdgeColor',[0,0,0]);
histogram(PCTr(sig1),'BinEdges',0:0.02:1,'Normalization','percentage',...
    'FaceAlpha',0.5,'FaceColor',[0,0,1],'EdgeAlpha',0.5,'EdgeColor',[0,0,0.7]);
histogram(PCTr(sig2),'BinEdges',1:0.02:2,'Normalization','percentage',...
    'FaceAlpha',0.5,'FaceColor',[1,0,0],'EdgeAlpha',0.5,'EdgeColor',[0.7,0,0]);

plot(PCT(:,1)*100,PCT(:,2)*100,'k.');
plot(PCT(sig1,1)*100,PCT(sig1,2)*100,'b.');
plot(PCT(sig2,1)*100,PCT(sig2,2)*100,'r.');
set(gca,'XLim',[15 65],'XTick',20:10:60,'YLim',[15 65],'YTick',20:10:60)
ylabel({'frequency decoding accuracy %','(pseudo-population)'},'FontName','Arial');
xlabel({'frequency decoding accuracy %','(session)'},'FontName','Arial');
title("Decoding performance per session",'FontName','Arial');
legend(["","No diff.","Pseudo > Sess","Sess > Pseudo"],'Box','off');
hold off;
[~,idx1] = min(allOS);
[~,idx2] = max(allOS);

subplot(2,3,3); hold on;
swarmchart(ones(nIter,1),squeeze(PCTall(:,idx1,1,1,1)).*100,'blue','o','SizeData',6);
swarmchart(ones(nIter,1).*2,squeeze(PCTall(:,idx1,1,1,2)).*100,'blue','o','SizeData',6);
set(gca,'XLim',[0 3],'XTick',[1 2],'XTickLabels',["Sess","Pseudo"],'YLim',[0 100],...
    'YTick',[0 50 100]);
ylabel('accuracy %','FontName','Arial');
title("Example session",'FontName','Arial');
hold off;

subplot(2,3,6); hold on;
swarmchart(ones(nIter,1),squeeze(PCTall(:,idx2,1,1,1)).*100,'red','o','SizeData',6);
swarmchart(ones(nIter,1).*2,squeeze(PCTall(:,idx2,1,1,2)).*100,'red','o','SizeData',6);
set(gca,'XLim',[0 3],'XTick',[1 2],'XTickLabels',["Sess","Pseudo"],'YLim',[0 100],...
    'YTick',[0 50 100]);
ylabel('accuracy %','FontName','Arial');
title("Example session",'FontName','Arial');
hold off;

% % Histogram splitting groups
% figure('Color',[1 1 1]); hold on;
% PCTr = PCT(:,1) ./ PCT(:,2);
% histogram(PCTr(~sig1 & ~sig2),'BinEdges',0:0.02:1,...
% 'FaceAlpha',0.5,'FaceColor',[.2,.2,.2],'EdgeAlpha',0.5,'EdgeColor',[0,0,0]);
% histogram(PCTr(sig1),'BinEdges',0:0.02:1,...
% 'FaceAlpha',0.5,'FaceColor',[0,0,1],'EdgeAlpha',0.5,'EdgeColor',[0,0,0.7]);
% histogram(PCTr(sig2),'BinEdges',1:0.02:2,...
% 'FaceAlpha',0.5,'FaceColor',[1,0,0],'EdgeAlpha',0.5,'EdgeColor',[0.7,0,0]);
% ylabel({'Number of sessions','(pseudo-population)'},'FontName','Arial');
% xlabel({"Ratio of frequency decoding accuracy","(Session : Pseudo-pop.)"},'FontName','Arial');
% title("Decoding performance per session",'FontName','Arial');
% set(gca,'XLim',[0.5 2]);

% % Which mice provide good per-session decoding?
% figure('Color',[1 1 1]);
% subplot(1,2,1); hold on;
% locSessNames = char(sessIDs(sig2));
% locSessNames = string(locSessNames(:,1:4));
% histogram(categorical(locSessNames))
% ylabel('number of sessions','FontName','Arial');
% xlabel('mouse number','FontName','Arial');
% title("Session source (sess>pseudo)",'FontName','Arial');
% hold off;
% subplot(1,2,2); hold on;
% locSessNames = char(sessIDs(sig1));
% locSessNames = string(locSessNames(:,1:4));
% histogram(categorical(locSessNames))
% ylabel('number of sessions','FontName','Arial');
% xlabel('mouse number','FontName','Arial');
% title("Session source (sess<pseudo)",'FontName','Arial');


% % What cell types are showing up in good decoding sessions
% figure('Color',[1 1 1]);
% hold on;
% for ns = 1 : length(sessIDs)
%     numNeu = sum(is.sessID(:,ns));
%     numPV = sum(is.sessID(:,ns) & is.PV);
%     numSOM = sum(is.sessID(:,ns) & is.SOM);
%     locAcc = squeeze(mean(PCTall(:,ns,1,1,1),1))*100;
%     if numPV > 0 && numSOM == 0
%         continue
%         % locRatio = numPV / numNeu;
%         % plot(locRatio,locAcc,'o','Color',colors.PV);
%         plot(numPV,locAcc,'o','Color',colors.PV);
%         % plot3(numPV,numNeu-numPV,locAcc,'o','Color',colors.PV);
%     elseif numSOM > 0 && numPV == 0
%         % locRatio = numSOM / numNeu;
%         % plot(locRatio,locAcc,'o','Color',colors.SOM);
%         plot(numSOM,locAcc,'o','Color',colors.SOM);
%         % plot3(numSOM,numNeu-numSOM,locAcc,'o','Color',colors.SOM);
%     end
% end
% xlabel('Inh:Exc Ratio','FontName','Arial');
% ylabel('Frequency decoding accuracy','FontName','Arial');
% % ylabel('number of EXC neurons','FontName','Arial');
% % xlabel('number of INH neurons','FontName','Arial');
% % zlabel('model accuracy','FontName','Arial');
% title("Session source (sess<pseudo)",'FontName','Arial');


%% DECODING ACCURACY BY NUMBER OF FREQUENCY SELECTIVE NEURONS PER SESSION

neuPerSess = zeros(size(is.sessID,2),1);
for ns = 1:length(sessIDs)
    neuPerSess(ns) = sum(is.sessID(:,ns) & is.allMod);
end
nonZeroSess = PCT(:,1) ~= 0;
locAcc = PCT(nonZeroSess,1);
locNums = neuPerSess(nonZeroSess);
% figure('Color',[1 1 1]); hold on;
% plot(locNums,locAcc,'ko','MarkerSize',4);
% mdl1 = fitlm(locNums,locAcc);
% if mdl1.ModelFitVsNullModel.Pvalue < pValThresh
%     % significant fit
%     xGrid = linspace(min(locNums), max(locNums), 100)';  
%     [y,ci] = predict(mdl1, xGrid);
%     plot(xGrid, y, 'r-', xGrid, ci(:,1), 'r:', xGrid, ci(:,2), 'r:', 'LineWidth', 1);
% end
% %set(gca,'XTick',1:6,'XTickLabels',0:0.2:1);
% ylabel('Frequency decoding accuracy (%)','FontName','Arial');
% xlabel('Number of frequency-selective units','FontName','Arial');
% title("Model performance within sessions",'FontName','Arial');
% text(100,0.3,strcat("p-val = ",num2str(mdl1.ModelFitVsNullModel.Pvalue,3)),'FontName','Arial');
% text(100,0.28,strcat("r-sq. = ",num2str(mdl1.Rsquared.Ordinary,3)),'FontName','Arial');

% AS SLIDING WINDOW AVERAGE
binWidth = 14; % number of neurons
binStarts = min(locNums):1:max(locNums)-binWidth; % starting positions (left hand bin edge)
allCIs = zeros(2,length(binStarts)); % all confidence intervals for each bin
allMeansFS = zeros(length(binStarts),1); % all mean accuracy for each bin
for nb = 1 : length(binStarts)
    locEdges = [binStarts(nb) binStarts(nb)+binWidth];
    locSessions = locNums >= locEdges(1) & locNums <= locEdges(2);
    allMeansFS(nb) = mean(locAcc(locSessions))*100;
    allCIs(:,nb) = mkCI(locAcc(locSessions))*100;
end
binAvgsFS = binStarts + (binWidth/2); % average value, center of bin
figure('Color',[1 1 1]); hold on;
patchY = cat(1,allCIs(2,:)',flip(allCIs(1,:))');
patchX = cat(1,binAvgsFS',flip(binAvgsFS)');
fill(patchX,patchY,[0.4 0 0.5],'FaceAlpha',0.3,'EdgeColor','none');
plot(binAvgsFS,allMeansFS,'-','Color',[0.4 0 0.5]);

% NOW COMPUTE USING ALL NEURONS AVAILABLE
neuPerSess = zeros(size(is.sessID,2),1);
for ns = 1:length(sessIDs)
    neuPerSess(ns) = sum(is.sessID(:,ns));
end
nonZeroSess = PCT(:,1) ~= 0;
locAcc = PCT(nonZeroSess,1);
locNums = neuPerSess(nonZeroSess);
binWidth = 20; % number of neurons
binStarts = min(locNums):1:max(locNums)-binWidth; % starting positions (left hand bin edge)
allCIs = zeros(2,length(binStarts)); % all confidence intervals for each bin
allMeansAA = zeros(length(binStarts),1); % all mean accuracy for each bin
for nb = 1 : length(binStarts)
    locEdges = [binStarts(nb) binStarts(nb)+binWidth];
    locSessions = locNums >= locEdges(1) & locNums <= locEdges(2);
    allMeansAA(nb) = mean(locAcc(locSessions))*100;
    allCIs(:,nb) = mkCI(locAcc(locSessions))*100;
end
binAvgsAA = binStarts + (binWidth/2); % average value, center of bin
patchY = cat(1,allCIs(2,:)',flip(allCIs(1,:))');
patchX = cat(1,binAvgsAA',flip(binAvgsAA)');
fill(patchX,patchY,[0.5 0 0.4],'FaceAlpha',0.3,'EdgeColor','none');
plot(binAvgsAA,allMeansAA,'--','Color',[0.5 0 0.4]);

% PLOT ALTERATIONS
set(gca,'XLim',[18.5 215],'YLim',[28.7 83.4]);
ylabel('Frequency decoding accuracy (%)','FontName','Arial');
xlabel('Number of neurons per session','FontName','Arial');
title("Model performance within sessions",'FontName','Arial');
legend({'','Frequency-Selective','','All Available'},'Box','off');

% SMALLER VERSION: BAR PLOT SHOWING DIFFERENCE FS > AA
[mn,mx] = bounds(intersect(binAvgsFS, binAvgsAA)); % min and max common values (x values)
idxFS = binAvgsFS >= mn & binAvgsFS <= mx; % index of common values in FS group
idxAA = binAvgsAA >= mn & binAvgsAA <= mx; % index of common values in AA group
locDiff = allMeansFS(idxFS) - allMeansAA(idxAA);
locMean = mean(locDiff);
locCI = mkCI(locDiff);
figure('Color',[1 1 1]); hold on;
bar(1,locMean,'FaceColor',[0.5 0 0.5]);
errorbar(1,locMean,locMean-locCI(1),'vertical','Color',[0 0 0]);
set(gca,'XTick',[],'XLim',[0 2]);

% TWO BAR PLOTS, FS vs AA
figure('Color',[1 1 1]); hold on;
bar(1,mean(allMeansFS(idxFS)),'FaceColor',[0.4 0 0.5]);
locCI = mkCI(allMeansFS(idxFS));
errorbar(1,mean(allMeansFS(idxFS)),mean(allMeansFS(idxFS))-locCI(1),'Color',[0 0 0]);
bar(2,mean(allMeansAA(idxAA)),'FaceColor',[0.6 0 0.4]);
locCI = mkCI(allMeansAA(idxAA));
errorbar(2,mean(allMeansAA(idxAA)),mean(allMeansAA(idxAA))-locCI(1),'Color',[0 0 0]);
set(gca,'XTick',[],'XLim',[0.2 2.8],'YLim',[30 60]);
ylabel('Frequency decoding accuracy (%)','FontName','Arial');
title("Model performance within sessions",'FontName','Arial');
legend({'Frequency-Selective','','All Available',''},'Box','off');





%% 

%% PLOT INSET FIGURE SHOWING TOP 3 PCs OF EACH CELL TYPE (DATASET 1)

targetSess = [42,59,83];
pcVisNum = 3; % number of PCs to visualize
figure('theme','light','color',[1 1 1]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
spx = [-10 0 10];
ntx = [1 2 3; 4 5 6; 7 8 9];

for ns = 1 : length(targetSess) % loop through example sessions
    for nr = 1 : 3 % loop through 3 repetitions train/test
        targetRep = nr;
        mdl = MDLall{targetRep,targetSess(ns),1,1,1}; 
        PCresps = zeros(trainTrial,5,pcVisNum);
        for np = 1 : pcVisNum
            for nf = 1 : nFreq
                PCresps(:,nf,np) = mdl.X(mdl.Y == freqVal(nf),np);
                % yMean = mean(PCresps(:,nf,np));
            end
        end
        nexttile(ntx(ns,nr)); hold on;
        for np = 1 : pcVisNum
            xVals = freqVal + spx(np); % frequency + xspacing  
            locResps = squeeze(PCresps(:,:,np));
            for nf = 1 : nFreq
                plot([xVals(nf) xVals(nf)],mkCI(locResps(:,nf)),'-','Color',[0 0 0]); % CI bar
                plot(xVals(nf),mean(locResps(:,nf)),'.','Color',[0 0 0]) % mean point
            end 
            switch np
                case 1
                    plot(xVals,mean(locResps,1),'-','Color',[0 0 0]);
                case 2
                    plot(xVals,mean(locResps,1),'--','Color',[0 0 0]);
                case 3
                    plot(xVals,mean(locResps,1),':','Color',[0 0 0]);
            end
        end
        % set(gca,'YLim',[-10 8],'YTick',-8:4:8,'XTick',freqVal,'XLim',[280 1120],'FontSize',15);
        set(gca,'YLim',[-5 5],'YTick',-4:2:4,'XTick',freqVal,'XLim',[280 1120],'FontSize',15);
        if nr == 1 
            ylabel("\DeltaF/F",'FontName','Arial');
        end
        if ns ~= 3
            set(gca,'XTickLabel',[]);
        end
        if ns == 3 && nr == 2
            xlabel("Frequency (Hz)",'FontName','Arial');
        end
        hold off
    end
end


% LEGEND IN SPOT 3
figure; hold on;
for nc = 1 : 3
plot(0,0,'.','Color',colors.(cellType(nc)),'MarkerSize',12);
end
plot([0 0],[0 0],'k-'); 
plot([0 0],[0 0],'k--');
plot([0 0],[0 0],'k:');
set(gca,'XLim',[1 2],'XColor','none','YColor','none','FontSize',12);
legend({"PV","SOM","EXC","PC 1","PC 2","PC 3"},'NumColumns',2,'box','off','FontSize',15);
% legend({"PV","SOM","EXC","PC 2","PC 3","PC 4"},'NumColumns',2,'box','off','FontSize',15);



%%  EXC/INH LOADING INTO EACH PC (DATASET 1)
figure('theme','light','color',[1 1 1]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
ntx = [1 2 3; 4 5 6; 7 8 9];

for ns = 1 : length(targetSess) % loop through example sessions
    for nr = 1 : 3 % loop through 3 repetitions train/test
        ctr1 = 1; %counter
        nexttile(ntx(ns,nr)); hold on;
        targetRep = nr;
        coefs = PCFall{targetRep,targetSess(ns),1,1,1};
        for np = 1 : 3
            locCoefs = coefs(:,np);
            locNeuID = Tc.identity(is.sessID(:,targetSess(ns)));
            locInh = locNeuID == "PV" | locNeuID == "SOM";
            locExc = locNeuID == "EXC";
            if any(locNeuID == "PV")
                swarmchart(ones(sum(locInh),1)*ctr1,locCoefs(locInh),'filled','o','MarkerFaceColor',...
                    colors.PV,'MarkerEdgeAlpha',0,'SizeData',8);
            elseif any(locNeuID == "SOM")
                swarmchart(ones(sum(locInh),1)*ctr1,locCoefs(locInh),'filled','o','MarkerFaceColor',...
                    colors.SOM,'MarkerEdgeAlpha',0,'SizeData',8);
            end
            ctr1 = ctr1 + 1; % inc counter
            swarmchart(ones(sum(locExc),1)*ctr1,locCoefs(locExc),'filled','o','MarkerFaceColor',...
                    colors.EXC,'MarkerEdgeAlpha',0,'SizeData',8);
            ctr1 = ctr1 + 1; % inc counter
        end
        set(gca,'XTick',[1.5 3.5 5.5],'XTickLabel',[])
        if ns == 3
            set(gca,'XTickLabels',{"PC1","PC2","PC3"});
        end
        if nr == 1 && ns == 2
            ylabel("Neuron PCA coeff. val.",'FontName','Arial');
        end
        if ns == 1
            title(strcat("Train Rep. #",num2str(nr)),'FontName','Arial');
        end
    end
end

%% COMPUTE: Avg neuron pairwise correlation in each session 
% DATASET 1
% (WARNING: ~1h20m run time)

disp("starting correlation analysis...");
allCorrVals = cell(length(sessIDs),1);
allPairNums = cell(length(sessIDs),1);
for ns = 1 : size(is.sessID,2) % loop through sessions
    tic;
    locNeuList = find(is.sessID(:,ns)); % local neuron list
    locnn = length(locNeuList);
    nPairs = locnn*(locnn-1)/2; % number of unique pairs from n neurons
    corrVals = zeros(nPairs,1);
    pairNums = zeros(nPairs,2);
    corrIdx = 1;
    disp(sessIDs(ns));
    % loop to make pairwise comparisons between all neurons
    for n1 = 1 : locnn-1 
        for n2 = n1 + 1 : locnn
            neu1 = locNeuList(n1); % neuron 1
            neu2 = locNeuList(n2); % neuron 2
            pairNums(corrIdx,:) = [neu1 neu2]; 
            % corrVals(corrIdx) = corr(Tc.responses{neu1},Tc.responses{neu2},'Type','Spearman');
            corrVals(corrIdx) = corr(Tc.lateResps{neu1},Tc.lateResps{neu2},'Type','Spearman');
            corrIdx = corrIdx + 1;
        end
    end
    allCorrVals{ns} = corrVals;
    allPairNums{ns} = pairNums;
    toc;
end
disp("finished correlation analysis!");

% % SAVE OUT
% cd('C:\Users\skich\Desktop\WORK');
save("dataset1_neuCorrValsPerSess_lateRespWindows.mat",'allCorrVals','allPairNums','-v7.3');

%% COMPUTE: Avg neuron pairwise correlation in each session 
% DATASET 2
% (WARNING: ~1h20m run time)

disp("starting correlation analysis...");
allCorrVals = cell(length(sessIDs),1);
allPairNums = cell(length(sessIDs),1);
for ns = 1 : size(is.sessID,2) % loop through sessions
    tic;
    locNeuList = find(is.sessID(:,ns)); % local neuron list
    locnn = length(locNeuList);
    nPairs = locnn*(locnn-1)/2; % number of unique pairs from n neurons
    corrVals = zeros(nPairs,1);
    pairNums = zeros(nPairs,2);
    corrIdx = 1;
    disp(sessIDs(ns));
    % loop to make pairwise comparisons between all neurons
    for n1 = 1 : locnn-1 
        for n2 = n1 + 1 : locnn
            neu1 = locNeuList(n1); % neuron 1
            neu2 = locNeuList(n2); % neuron 2
            idx1 = 1;
            idx2 = 1;
            numTrials = min([max(size(Tc.responses{neu1})) max(size(Tc.responses{neu2}))]);
            locResps = zeros(numTrials,2);
            respIdx = 1;
            for na = 1 : 5 % loop amplitude
                for nf = 1 : 5 % loop frequency
                    log1 = Tc.outlierTrials{neu1}{na,nf}; % logical array of outlier trials (1 = outlier)
                    log2 = Tc.outlierTrials{neu2}{na,nf};
                    for nn = 1 : length(log1) % loop through trials
                        if log1(nn) == 0 && log2(nn) == 0
                            % both trials are non-outliers
                            % locResps = cat(1,locResps,[Tc.responses{neu1}(idx1) Tc.responses{neu2}(idx2)]);
                            % locResps(respIdx,:) = [Tc.responses{neu1}(idx1) Tc.responses{neu2}(idx2)]; % 0-1sec dF/F avg.s
                            locResps(respIdx,:) = [Tc.lateResps{neu1}(idx1) Tc.lateResps{neu2}(idx2)];
                            respIdx = respIdx + 1;
                            idx1 = idx1 + 1;
                            idx2 = idx2 + 1;
                        elseif log1(nn) == 1 && log2(nn) == 0
                            % neu1 trial is outlier
                            % that means... neu2 was assigned a value in responses, which does not correspnd to neu1
                            % therefore, we should increase the index for neu2, so we skip the response value where
                            % there is no corresponding value in neuron 1
                            idx2 = idx2 + 1;
                        elseif log1(nn) == 0 && log2(nn) == 1
                            % neu2 trial is outlier, same logic as above
                            idx1 = idx1 + 1;
                        elseif log1(nn) == 1 && log2(nn) == 1
                            % neu1 AND neu2 trials are outliers
                        end
                    end
                end
            end
            pairNums(corrIdx,:) = [neu1 neu2]; 
            corrVals(corrIdx) = corr(locResps(:,1),locResps(:,2),'Type','Spearman');
            corrIdx = corrIdx + 1;
        end
    end
    allCorrVals{ns} = corrVals;
    allPairNums{ns} = pairNums;
    toc;
end
disp("finished correlation analysis!");

% % SAVE OUT
% cd('C:\Users\skich\Desktop\WORK');
% save("dataset2_neuCorrValsPerSess_lateRespWindows.mat",'allCorrVals','allPairNums','-v7.3');

%% Alternate analysis

% FAS_neuCorrByFreqDecoding.m

% -------------------------------------------------------------------------

%% plots for Sess Decoding vs. neuron correlations

% cd('C:\Users\skich\Desktop\WORK');
% load("dataset2_neuCorrValsPerSess.mat"); % 'allCorrVals' & 'allPairNums'

% Filter out trained sessions
trainedSess = false(size(is.sessID,2),1);
for ns = 1 : size(is.sessID,2)
    locNeuList = find(is.sessID(:,ns));
    trainedSess(ns) = Tc.isTrained(locNeuList(1));
end

% PLOT ALL SESS, CORR BY ACCURACY
figure('Color',[1 1 1]); hold on;
PCTloc = squeeze(mean(wPCA.PCTall(:,~trainedSess,1,1,1),1))' .* 100;
meanCorr = cellfun(@mean,allCorrVals(~trainedSess));
plot(meanCorr,PCTloc,'k.');
ylabel('Frequency decoding accuracy %','FontName','Arial');
xlabel('Mean neuron corr value','FontName','Arial');
title({"Neuron correlation by frequency decoding accuracy","(per session)"},'FontName','Arial');
mdl = fitlm(meanCorr,PCTloc);
xGrid = linspace(min(meanCorr), max(meanCorr), 100)';  
[y,ci] = predict(mdl,xGrid);
plot(xGrid, y, 'r-', xGrid, ci(:,1), 'r:', xGrid, ci(:,2), 'r:', 'LineWidth', 1);
text(.08,20,strcat("p-val = ",num2str(mdl.ModelFitVsNullModel.Pvalue,3)),'FontName','Arial');
text(.08,18,strcat("r-sq. = ",num2str(mdl.Rsquared.Ordinary,3)),'FontName','Arial');
hold off;


% PLOT sess corr val by accuracy (ONLY SELECTIVE NEURONS)
figure('Color',[1 1 1]); hold on;
locMeanCorr = zeros(size(allPairNums,1),1);
for ns = 1 : size(allPairNums,1) % loop sessions
    locIdx = all(is.allMod(allPairNums{ns}),2); % neuron pairs where both are selective
    locVals = allCorrVals{ns}(locIdx);
    locMeanCorr(ns) = mean(locVals); % corr values from selective pairs
end
locMeanCorr = locMeanCorr(~trainedSess);
plot(locMeanCorr,PCTloc,'k.');
ylabel('Frequency decoding accuracy %','FontName','Arial');
xlabel('Mean neuron corr value','FontName','Arial');
title({"Neuron correlation by frequency decoding accuracy","freq. sel. neu. pairs only"},'FontName','Arial');
mdl = fitlm(locMeanCorr,PCTloc);
xGrid = linspace(min(locMeanCorr), max(locMeanCorr), 100)';  
[y,ci] = predict(mdl,xGrid);
plot(xGrid, y, 'r-', xGrid, ci(:,1), 'r:', xGrid, ci(:,2), 'r:', 'LineWidth', 1);
text(.2,24,strcat("p-val = ",num2str(mdl.ModelFitVsNullModel.Pvalue,3)),'FontName','Arial');
text(.2,22,strcat("r-sq. = ",num2str(mdl.Rsquared.Ordinary,3)),'FontName','Arial');
%set(gca,'XLim',[-0.05 0.2],'YLim',[15 70]);
hold off;

% split pairwise correlations up by cell type (PV|SOM|EXC)
allPairTypes = cell(size(allCorrVals));
for ns = 1 : size(allPairNums,1)
    allPairTypes{ns} = zeros(length(allPairNums{ns}),1);
    for np = 1 : length(allPairNums{ns})
        locPairID = Tc.identity(allPairNums{ns}(np,:));
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

% PLOT sess corr val (for each pair type) by accuracy

figure('Color',[1 1 1]);
tiledlayout(2,3);
pairNames = {"EXC:EXC", "EXC:PV", "EXC:SOM", "PV:PV", "SOM:SOM"};

% EXC:EXC
npt = 1;
locPairVals = NaN(size(allPairNums,1),1);
locMeanCorr = NaN(size(allPairNums,1),1);
locInhType = zeros(size(allPairNums,1),1);
for ns = 1 : size(allPairNums,1) % loop sessions
    locIdx = all(is.allMod(allPairNums{ns}),2); % neuron pairs where both are selective
    locVals = allCorrVals{ns}(allPairTypes{ns} == npt & locIdx);
    locPairVals(ns) = mean(locVals);
    locMeanCorr(ns) = mean(locVals); % corr values from selective pairs
    if any(allPairTypes{ns}==2) % PV present
        locInhType(ns) = 1;
    elseif any(allPairTypes{ns}==3) % SOM present
        locInhType(ns) = 2;
    end
end
locMeanCorr = locMeanCorr(~trainedSess);
abvThresh = locMeanCorr > 0.3;
locMeanCorr(abvThresh) = [];
PCTloc0 = PCTloc(~abvThresh);
locInhType(abvThresh) = [];
nexttile; hold on;
locCorrPV = []; locAccPV = [];
locCorrSOM = []; locAccSOM = [];
for np = 1 : length(locMeanCorr)
    if locInhType(np) == 1
        % plot(locMeanCorr(np),PCTloc0(np),'.','Color',colors.PV); % PV
        plot(locMeanCorr(np),PCTloc0(np),'.','Color',[0 0 0]); % PV
        locCorrPV = [locCorrPV; locMeanCorr(np)]; 
        locAccPV = [locAccPV; PCTloc0(np)];
    elseif locInhType(np) == 2
        % plot(locMeanCorr(np),PCTloc0(np),'.','Color',colors.SOM); % SOM
        plot(locMeanCorr(np),PCTloc0(np),'.','Color',[0 0 0]); % SOM
        locCorrSOM = [locCorrSOM; locMeanCorr(np)]; 
        locAccSOM = [locAccSOM; PCTloc0(np)];
    end
end
ylabel('Frequency decoding accuracy %','FontName','Arial');
xlabel('Mean neuron corr value','FontName','Arial');
title({"Neuron correlation by frequency decoding accuracy",pairNames{npt}},'FontName','Arial');

% FIT ALL POINTS
mdl = fitlm(locMeanCorr,PCTloc0);
xGrid = linspace(min(locMeanCorr), max(locMeanCorr), 100)'; 
[y,ci] = predict(mdl,xGrid);
plot(xGrid, y, 'r-', xGrid, ci(:,1), 'r:', xGrid, ci(:,2), 'r:', 'LineWidth', 1);
text(-.2,62.5,strcat("p-val = ",num2str(mdl.ModelFitVsNullModel.Pvalue,3)),'FontName','Arial');
text(-.2,57.5,strcat("r-sq. = ",num2str(mdl.Rsquared.Ordinary,3)),'FontName','Arial');

% % FIT PV ponts
% mdl = fitlm(locCorrPV,locAccPV);
% xGrid = linspace(min(locCorrPV), max(locCorrPV), 100)'; 
% [y,ci] = predict(mdl,xGrid);
% plot(xGrid, y, '-', xGrid, ci(:,1), ':', xGrid, ci(:,2), ':', 'LineWidth', 1,'Color',colors.PV);
% text(-.2,62.5,strcat("p-val = ",num2str(mdl.ModelFitVsNullModel.Pvalue,3)),'FontName','Arial','Color',colors.PV);
% text(-.2,57.5,strcat("r-sq. = ",num2str(mdl.Rsquared.Ordinary,3)),'FontName','Arial','Color',colors.PV);
% % FIT SOM ponts
% mdl = fitlm(locCorrSOM,locAccSOM);
% xGrid = linspace(min(locCorrSOM), max(locCorrSOM), 100)'; 
% [y,ci] = predict(mdl,xGrid);
% plot(xGrid, y, '-', xGrid, ci(:,1), ':', xGrid, ci(:,2), ':', 'LineWidth', 1,'Color',colors.SOM);
% text(-.2,52.5,strcat("p-val = ",num2str(mdl.ModelFitVsNullModel.Pvalue,3)),'FontName','Arial','Color',colors.SOM);
% text(-.2,47.5,strcat("r-sq. = ",num2str(mdl.Rsquared.Ordinary,3)),'FontName','Arial','Color',colors.SOM);

set(gca,'XLim',[-0.22 0.3],'YLim',[15 65]);
hold off;

% OTHER PAIRWISE TYPES
for npt = 2 : 5 % num pair types (EXC:EXC, EXC:PV, EXC:SOM, PV:PV, SOM:SOM)
    locPairVals = NaN(size(allPairNums,1),1);
    locMeanCorr = NaN(size(allPairNums,1),1);
    for ns = 1 : size(allPairNums,1) % loop sessions
        locIdx = all(is.allMod(allPairNums{ns}),2); % neuron pairs where both are selective
        locVals = allCorrVals{ns}(allPairTypes{ns} == npt & locIdx);
        locPairVals(ns) = mean(locVals);
        locMeanCorr(ns) = mean(locVals); % corr values from selective pairs
    end
    locMeanCorr = locMeanCorr(~trainedSess);
    nexttile; hold on;
    % plot(locPairVals,PCTloc,'k.');
    plot(locMeanCorr,PCTloc,'k.');
    ylabel('Frequency decoding accuracy %','FontName','Arial');
    xlabel('Mean neuron corr value','FontName','Arial');
    title({"Neuron correlation by frequency decoding accuracy",pairNames{npt}},'FontName','Arial');
    % mdl = fitlm(locPairVals,PCTloc);
    mdl = fitlm(locMeanCorr,PCTloc);
    % xGrid = linspace(min(locPairVals), max(locPairVals), 100)';  
    xGrid = linspace(min(locMeanCorr), max(locMeanCorr), 100)'; 
    [y,ci] = predict(mdl,xGrid);
    plot(xGrid, y, 'r-', xGrid, ci(:,1), 'r:', xGrid, ci(:,2), 'r:', 'LineWidth', 1);
    text(-.2,62.5,strcat("p-val = ",num2str(mdl.ModelFitVsNullModel.Pvalue,3)),'FontName','Arial');
    text(-.2,57.5,strcat("r-sq. = ",num2str(mdl.Rsquared.Ordinary,3)),'FontName','Arial');
    set(gca,'XLim',[-0.22 0.3],'YLim',[15 65]);
    hold off;
end






%% PCA VS. no PCA

% save('modelResults_dataset2_perSession_PCAvNO.mat','wPCA','woPCA','-v7.3');
% load('modelResults_dataset2_perSession_PCAvNO.mat');


wPCA.PCT = squeeze(mean(wPCA.PCTall,1));
woPCA.PCT = squeeze(mean(woPCA.PCTall,1));

n1 = 2; % session [1], pseudo-pop [2]

% test for differences between wPCA and woPCA per session
allPVals = zeros(ns,1);
allOS = zeros(ns,1); % observed stats
for n = 1 : size(wPCA.PCTall,2)
    array1 = wPCA.PCTall(:,n,1,1,n1);
    array2 = woPCA.PCTall(:,n,1,1,n1);
    [pVal, obsStat, permStats] = permTest2sample(array1, array2, 5000);
    allPVals(n) = pVal;
    allOS(n) = obsStat;
end
sig1 = allPVals < pValThresh & allOS < 0;
sig2 = allPVals < pValThresh & allOS > 0;

% SCATTER (SESS PCA vs. noPCA) + 2 EXAMPLES
figure('Color',[1 1 1]);
subplot(2,3,[1 2 4 5]); hold on;
plot([0 100],[0 100],'k:');
plot(wPCA.PCT(:,n1)*100,woPCA.PCT(:,n1)*100,'k.');
plot(wPCA.PCT(sig1,n1)*100,woPCA.PCT(sig1,n1)*100,'b.');
plot(wPCA.PCT(sig2,n1)*100,woPCA.PCT(sig2,n1)*100,'r.');

%set(gca,'XLim',[15 65],'YLim',[25 45])
ylabel({'frequency decoding accuracy %','(without PCA)'},'FontName','Arial');
xlabel({'frequency decoding accuracy %','(with PCA)'},'FontName','Arial');
title("PCA effect on per session frequency decoding",'FontName','Arial');
legend(["","No diff.","noPCA > PCA","PCA > noPCA"],'Box','off');
hold off;

subplot(2,3,3); hold on;
[~,idx1] = min(allOS); % woPCA >> PCA
swarmchart(ones(nIter,1),squeeze(wPCA.PCTall(:,idx1,1,1,n1)).*100,'blue','o','SizeData',6);
swarmchart(ones(nIter,1).*2,squeeze(woPCA.PCTall(:,idx1,1,1,n1)).*100,'blue','o','SizeData',6);
set(gca,'XLim',[0 3],'XTick',[1 2],'XTickLabels',["PCA","no PCA"],'YLim',[0 100],...
    'YTick',[0 50 100]);
ylabel('accuracy %','FontName','Arial');
title("Example session",'FontName','Arial');
hold off;

subplot(2,3,6); hold on;
[~,idx2] = max(allOS); % PCA >> woPCA
swarmchart(ones(nIter,1),squeeze(wPCA.PCTall(:,idx2,1,1,n1)).*100,'red','o','SizeData',6);
swarmchart(ones(nIter,1).*2,squeeze(woPCA.PCTall(:,idx2,1,1,n1)).*100,'red','o','SizeData',6);
set(gca,'XLim',[0 3],'XTick',[1 2],'XTickLabels',["PCA","no PCA"],'YLim',[0 100],...
    'YTick',[0 50 100]);
ylabel('accuracy %','FontName','Arial');
title("Example session",'FontName','Arial');
hold off;



%% HOW DOES AMPLITUDE EFFECT PERFORMANCE

load('modelAccuracy_dataset2_perSession_allAmps.mat');

figure;
swarmchart(repmat(ampValAll',94,1),squeeze(PCT_byAmp(:,:,1))');



%% PLOTS
 
 
% % Decoding accuracy by avg sel str
% locSelStr = zeros(length(sessIDs),1);
% for ns = 1 : length(sessIDs)
%     locSelStr(ns) = mean(Tc.selectivityMetric(is.sessID(:,ns)));
% end
% locAcc = mean(squeeze(PCTall(:,:,6)),1);
% figure('Color',[1 1 1]); hold on;
% plot(locSelStr,locAcc,'ko','MarkerSize',4);
% 
% mdl2 = fitlm(locSelStr,locAcc);
% if mdl2.ModelFitVsNullModel.Pvalue < pValThresh
%     % significant fit
%     xGrid = linspace(min(locSelStr), max(locSelStr), 100)';  
%     [y,ci] = predict(mdl2, xGrid);
%     plot(xGrid, y, 'r-', xGrid, ci(:,1), 'r:', xGrid, ci(:,2), 'r:', 'LineWidth', 1);
% end
% ylabel('Frequency decoding accuracy (%)','FontName','Arial');
% xlabel('Avg. Sel. Str. per session','FontName','Arial');
% title("Model performance by Freq. Sel. Str.",'FontName','Arial');
% text(.42,.14,strcat("p-val = ",num2str(mdl2.ModelFitVsNullModel.Pvalue,3)),'FontName','Arial');
% text(.42,.08,strcat("r-sq. = ",num2str(mdl2.Rsquared.Ordinary,3)),'FontName','Arial');
 
 
% % model performance scatter (gamma = 0, x-axis) gamma=1 y axis
% avgPerfG0 = mean(squeeze(PCTall(:,:,1)),1);
% avgPerfG1 = mean(squeeze(PCTall(:,:,6)),1);
% figure('Color',[1 1 1]); hold on;
% plot([0 1],[0 1],'r--');
% plot(avgPerfG0,avgPerfG1,'ko','MarkerSize',4);
% set(gca,'XLim',[0.08 0.7],'YLim',[0.08 0.7]);
% ylabel('Frequency decoding accuracy (gamma = 1)','FontName','Arial');
% xlabel('Frequency decoding accuracy (gamma = 0)','FontName','Arial');
% title("Model performance per session",'FontName','Arial');
 
 
% % Decoding accuracy by neu num for 80% CVE
% neuPerSess = sum(is.sessID,1);
% locACVE = mean(squeeze(CVEall(:,:,6)),1);
% locAcc = mean(squeeze(PCTall(:,:,6)),1);
% figure('Color',[1 1 1]); hold on;
% plot(locACVE,locAcc,'ko','MarkerSize',4);
% 
% mdl3 = fitlm(locACVE,locAcc);
% if mdl3.ModelFitVsNullModel.Pvalue < pValThresh
%     % significant fit
%     xGrid = linspace(min(locACVE), max(locACVE), 100)';  
%     [y,ci] = predict(mdl3, xGrid);
%     plot(xGrid, y, 'r-', xGrid, ci(:,1), 'r:', xGrid, ci(:,2), 'r:', 'LineWidth', 1);
% end
% ylabel('Frequency decoding accuracy (%)','FontName','Arial');
% xlabel('Number of neurons per session to explain 80% of variance','FontName','Arial');
% title("Model performance by Freq. Sel. Str.",'FontName','Arial');
% text(7,.62,strcat("p-val = ",num2str(mdl3.ModelFitVsNullModel.Pvalue,3)),'FontName','Arial');
% text(7,.58,strcat("r-sq. = ",num2str(mdl3.Rsquared.Ordinary,3)),'FontName','Arial');


% % Model performance all sessions, gamma 1/0
% figure('Color',[1 1 1]); hold on;
% [locVal,locIdx] = sort(mean(squeeze(PCTall(:,:,6)),1),'descend'); 
% plot(locVal,'ro','MarkerSize',4);
% plot(sort(mean(squeeze(PCTall(:,:,1)),1),'descend'),'bo','MarkerSize',4);
% %set(gca,'XTick',1:6,'XTickLabels',0:0.2:1);
% ylabel('Frequency decoding accuracy (%)','FontName','Arial');
% xlabel('Session number','FontName','Arial');
% title("Model performance per session",'FontName','Arial');
% legend({'gamma = 1','gamma = 0'},'Box','off');



% % Model performance all sessions, gamma 1/0 matched session id
% figure('Color',[1 1 1]); hold on;
% [locVal,locIdx] = sort(mean(squeeze(PCTall(:,:,1)),1),'descend');
% plot(locVal,'bo','MarkerSize',4);
% plot(mean(squeeze(PCTall(:,locIdx,6)),1),'ro','MarkerSize',4);
% legend({'gamma = 1','gamma = 0'},'Box','off');
% ylabel('Frequency decoding accuracy (%)','FontName','Arial');
% xlabel('Session number','FontName','Arial');
% title("Model performance per session",'FontName','Arial');
% legend({'gamma = 1','gamma = 0'},'Box','off');



% % per by gamma
% figure('Color',[1 1 1]); hold on;
% xSP = [-0.1,0,0.1];
% for ns = 1 : length(cellType)
%     allCIs =  zeros(2,length(gammaSteps));
%     for ng = 1 : length(gammaSteps)
%         allCIs(:,ng) = mkCI(squeeze(PCTall(:,1,ns,ng,1)));
%         plot([ng ng]+xSP(ns), allCIs(:,ng)','-','Color',colors.(cellType(ns)));
%         plot(ng+xSP(ns),mean(squeeze(PCTall(:,1,ns,ng,1))),'.','Color',colors.(cellType(ns)));
%     end
% end
% set(gca,'XTick',1:6,'XTickLabels',0:0.2:1);
% ylabel('Frequency decoding accuracy (%)','FontName','Arial');
% xlabel('Gamma value','FontName','Arial');
% title("Model performance by gamma",'FontName','Arial');


% % per by delta
% figure('Color',[1 1 1]); hold on;
% xSP = [-0.1,0,0.1];
% for ns = 1 : length(cellType)
%     allCIs =  zeros(2,length(deltaSteps));
%     for ng = 1 : length(deltaSteps)
%         allCIs(:,ng) = mkCI(squeeze(PCTall(:,1,ns,1,ng)));
%         plot([ng ng]+xSP(ns), allCIs(:,ng)','-','Color',colors.(cellType(ns)));
%         plot(ng+xSP(ns),mean(squeeze(PCTall(:,1,ns,1,ng))),'.','Color',colors.(cellType(ns)));
%     end
% end
% set(gca,'XTick',1:6,'XTickLabels',0:0.2:1);
% ylabel('Frequency decoding accuracy (%)','FontName','Arial');
% xlabel('Delta value','FontName','Arial');
% title("Model performance by gamma",'FontName','Arial');



%% Cellular makeup of good decoding sessions

% mouseNumPV = ["m172","m270","m627","m630","m669","m748","m792","m806","m810"];
% mouseNumSOM = ["m126","m373","m418","m422","m468","m469","m699"];





























