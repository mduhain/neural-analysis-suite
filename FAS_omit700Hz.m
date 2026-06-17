% FAS_omit700Hz
%
%

%% load

% CHANGE DIRECTORIES
% cd('C:\Users\ramirezlab\Box\Tactile_Synchrony\2P-Data\Cleaned'); % LAB MACHINES
% cd('Y:\Michael\Haptics-2P\Reports'); % LAB MACHINES
cd('C:\Users\skich\Box\Tactile_Synchrony\2P-Data\Cleaned');

tic; load('dataset1_lightweight.mat'); toc;
%tic; load("dataset2_lightweight.mat"); toc;
tic; load("dataset2_lightweight_final.mat"); toc;

% CELL TYPE SPECIFIC COLORS RGB-TRIPLETS
colors = struct();
colors.PV = [0.2 0.62 0.84];
colors.SOM = [0.9 0.55 0.1];
colors.Green = [0.2 0.2 0.2];
colors.EXC = [0.2 0.2 0.2];
colors.shift = 0.8;
colors.Resp = [0 0.4470 0.7410];
colors.Sel = [0.4660 0.6740 0.1880];
colors.Amp = [1,0,0; 0.5,0,0; 0.3,0.3,0.3; 0,0.5,0.5; 0,1,1]; %RED TO CYAN


%% DECLARE ACTIVE DATASET


Tc = T1; % <----

clear T1 T2 % remove others from memory to save RAM
cellType = unique(Tc.identity);
nCellType = length(cellType);
freqValAll = unique(Tc.frequencies{1});

% CHECK FOR MULTIPLE AMPLITUDES
if any(strcmp(fieldnames(Tc),'amplitudes'))
    multiAmp = true;
    ampValAll = unique(Tc.amplitudes{1});
    nAmp = length(ampValAll);
    freqVal = freqValAll(1:end);
    nFreq = length(freqVal);
else
    multiAmp = false;
    freqVal = freqValAll([2 3 5 6]); % skip 100 & 700 HZ
    nFreq = length(freqVal);
end

disp("Dataset Initiated");


%% filter data 

minTrialNum = 25; % 25 (dataset1) | 15 (dataset2)

minTrialCountsF = zeros(size(Tc,1),1); %for frequency alone
minTrialCounts = zeros(size(Tc,1),1); %for freq X amp

for i = 1:size(Tc,1)
    % Extract values
    cellFreq = Tc.frequencies{i};
    cellResp = Tc.responses{i};
    if multiAmp
        cellAmp = Tc.amplitudes{i};
    end

    %  Map each value in cellFreq to its index in freqVal
    [~, locF] = ismember(cellFreq,freqValAll);
    if multiAmp
        [~, locA] = ismember(cellAmp,ampValAll);
    end

    % Count occurrences of trials at each frequency
    if multiAmp %DATASET2
        jointCounts = accumarray([locF,locA],1,[length(freqValAll),length(ampValAll)]);
        minTrialCounts(i) = min(jointCounts,[],'all');
        if minTrialCounts(i) < minTrialNum
            % too few trials, exclude neuron
            Tc.frequencies{i} = [];
            Tc.responses{i} = [];
            Tc.amplitudes{i} = [];
        end

    else % DATASET1
        trialCountsF = accumarray(locF,1,[length(freqValAll),1]);
        % minimum trial number requirement for 300-1100 Hz
        minTrialCountsF(i) = min(trialCountsF(2:end));
        if  minTrialCountsF(i) >= minTrialNum
            % % exclude 100 Hz
            % Tc.frequencies{i} = cellFreq(locF>1); 
            % Tc.responses{i} = cellResp(locF>1);

            % exclude 100 Hz & 700 Hz
            Tc.frequencies{i} = cellFreq(locF~=1 & locF~=4); 
            Tc.responses{i} = cellResp(locF~=1 & locF~=4);
        else
            Tc.frequencies{i} = [];
            Tc.responses{i} = [];
        end
    end
    
end

disp('minTrialCounts')
if multiAmp
    tabulate(categorical(minTrialCounts))
else
    tabulate(categorical(minTrialCountsF))
end

% remove empty cells
Tc(cellfun(@isempty,Tc.responses),:) = [];
nCellTotal = size(Tc,1);

% BUILD IDENTIFICATION STRUCT 'IS'
is.PV = false(nCellTotal,1); % PV neuron
is.SOM = false(nCellTotal,1); % SOM neuron
is.EXC = false(nCellTotal,1); % EXC neuron
sessIDs = unique(Tc.sessID); % unique session IDs
mouseNums = unique(extractBefore(sessIDs, "_")); % unique mouse numbers
is.sessID = false(size(Tc,1),length(sessIDs)); % logical arrays: cells from each session(s)
is.mouseNum = false(size(Tc,1),length(mouseNums)); % logical arrays: cells from each mouse

for nc = 1 : nCellTotal % ASSIGN LOGICAL INDEXES
    is.sessID(nc,Tc.sessID(nc) == sessIDs) = true; % Session ID
    is.mouseNum(nc,extractBefore(Tc.sessID(nc),'_') == mouseNums) = true; % Mouse Number
    % Cell type
    switch string(Tc.identity(nc))
        case "EXC"
            is.EXC(nc) = true;
        case "PV"
            is.PV(nc) = true;
        case "SOM"
            is.SOM(nc) = true;
    end
end

disp('cell identity')
tabulate(Tc.identity)

%% calculate selectivity metric: anova effect size
% Pre-allocate columns in the table 

pValThresh = 0.05; % P VALUE THRESHOLD

Tc.responsivityEffectSize = zeros(nCellTotal, 1);
Tc.responsivityPVal = zeros(nCellTotal, 1);

if multiAmp == true
    Tc.anovaFreqP = zeros(nCellTotal,1);
    Tc.anovaAmpP = zeros(nCellTotal,1);
    Tc.anovaMixP = zeros(nCellTotal,1);
    % ADD STRUCT FIELDS TO 'IS'
    is.FreqMod = false(nCellTotal,1); % frequency modulated ONLY
    is.AmpMod = false(nCellTotal,1); % amplitude modulated ONLY
    is.DualMod = false(nCellTotal,1); % modulated by frequency AND amplitude
    is.InterMod = false(nCellTotal,1); % interaction effect
elseif multiAmp == false
    Tc.selectivityEffectSize = zeros(nCellTotal, 1);
    Tc.selectivityPVal = zeros(nCellTotal, 1);
end

for i = 1:nCellTotal
    cellFreq = Tc.frequencies{i};
    cellResp = Tc.responses{i};
    if multiAmp == true
        cellAmp = Tc.amplitudes{i};
        respPval = zeros(nFreq,nAmp);
        respES = zeros(nFreq,nAmp);

        % responsiveness (ttest and effect size)
        for nf = 1 : nFreq
            for na = 1 : nAmp
                tmp = cellResp(cellFreq == freqVal(nf) & cellAmp == ampValAll(na));
                [~,respPval(nf,na)] = ttest(tmp);
                respES(nf,na) = mean(tmp)/std(tmp);
            end
        end
        % select most responsive cells
        [minPval,minPvalId] = min(respPval,[],'all');
        Tc.responsivityEffectSize(i) = respES(minPvalId);
        Tc.responsivityPVal(i) = minPval;

        % Custom anovan 2 way
        [pVal, tbl, ~] = anovan(cellResp, {cellFreq, cellAmp}, 'model',...
            'interaction','varnames',{'FreqFactor','AmpFactor'},'display', 'off');
        Tc.anovanFreqP(i) = pVal(1); % Freq Effect
        Tc.anovanAmpP(i) = pVal(2); % Amp Effect
        Tc.anovanMixP(i) = pVal(3); % Interaction Effect
        % populate 'is' struct fields
        if  pVal(1) < pValThresh && pVal(2) > pValThresh && pVal(3) > pValThresh
            is.FreqMod(i) = true; % significant frequency effect ONLY, no interaction
        end
        if  pVal(1) > pValThresh && pVal(2) < pValThresh && pVal(3) > pValThresh
            is.AmpMod(i) = true; % significant amplitude effect ONLY, no interaction
        end
        if pVal(1) < pValThresh && pVal(2) < pValThresh  && pVal(3) > pValThresh
            is.DualMod(i) = true; % significant frequency AND amplitude effect, no interaction
        end
        if pVal(3) < pValThresh
            is.InterMod(i) = true; % significant interaction effect, any combination
        end
        
    elseif multiAmp == false
        % perform t-test for responsivity at each frequency
        respPval = zeros(nFreq,1);
        respES = zeros(nFreq,1);
        for j = 1:nFreq
            tmp = cellResp(cellFreq == freqVal(j));
            [~,respPval(j)] = ttest(tmp);
            respES(j) = mean(tmp)/std(tmp);
        end
        % select most responsive cells
        [minPval,minPvalId] = min(respPval);
        Tc.responsivityEffectSize(i) = respES(minPvalId);
        Tc.responsivityPVal(i) = minPval;
    
        % Perform one-way ANOVA for selectivity
        [p_val, tbl] = anova1(cellResp, cellFreq, 'off'); % 'off' suppresses the plot
        
        % Calculate eta-squared: SS_between / SS_total
        % tbl{2,2} is the Sum of Squares for the groups (Model)
        % tbl{4,2} is the Total Sum of Squares
        Tc.selectivityEffectSize(i) = tbl{2,2} / tbl{4,2};
        Tc.selectivityPVal(i) = p_val;
    end
end

is.responsive = Tc.responsivityPVal < pValThresh;
if multiAmp == false
    is.selective = Tc.selectivityPVal < pValThresh;
end

fprintf(newline);

% display count values
if multiAmp == true
    for n = 1 : length(cellType)
        loc1 = (is.(cellType(n)) & is.FreqMod);
        disp(strcat(cellType(n)," frequency effect: ",num2str(sum(loc1))," neurons"));
        loc2 = (is.(cellType(n)) & is.AmpMod);
        disp(strcat(cellType(n)," amplitude effect: ",num2str(sum(loc2))," neurons"));
        loc3 = (is.(cellType(n)) & is.DualMod);
        disp(strcat( cellType(n)," dual effect: ", num2str(sum(loc3))," neurons"));
        loc4 = (is.(cellType(n)) & is.InterMod);
        disp(strcat( cellType(n)," interaction effect: ", num2str(sum(loc4))," neurons"));
    end
    is.allMod = (is.FreqMod | is.DualMod | is.InterMod);
elseif multiAmp == false
    is.FreqMod = Tc.selectivityPVal < pValThresh;
    for n = 1 : length(cellType)
        loc1 = (is.(cellType(n)) & is.FreqMod);
        disp(strcat(cellType(n)," frequency effect: ",num2str(sum(loc1))," neurons"));
    end
end

fprintf(newline);
disp('Selectivity & ANOVAs calculated');

%% examine selectivity of each cell type 
figure('theme','light');

responsivityPercentage = zeros(nCellType,1);
selectivityPercentage = zeros(nCellType,1);
if multiAmp == true
    Tc.selectivityPVal = Tc.anovanFreqP;
end
for i = 1:length(cellType)
    idx = Tc.identity == cellType(i);    
    responsivityPercentage(i) = sum(Tc.responsivityPVal(idx)<0.05)/sum(idx)*100;
    selectivityPercentage(i) = sum(Tc.selectivityPVal(idx)<0.05)/sum(idx)*100;
end
nexttile;
bar(responsivityPercentage);
xticklabels(cellType); ylabel('Responsive cells %');ylim([0 100])
nexttile;
bar(selectivityPercentage);
xticklabels(cellType); ylabel('Selective cells %');ylim([0 100])

% examine cdf of each cell type
nexttile
pThr = 0.05;
for j = 1:length(cellType)
    idx = Tc.identity == cellType(j) & Tc.selectivityPVal < pThr;

    % calculate pdf of effect sizes for each cell type, 
    [pdfEffectSize, xValues] = ksdensity(Tc.selectivityEffectSize(idx),'Support',[0 1]);
    hold on;
    plot(xValues, pdfEffectSize, 'DisplayName', char(cellType(j)), 'LineWidth', 2);
end
xlabel('Selectivity Effect Size');
ylabel('Probability Density');
title('PDF of Effect Sizes by Cell Type');
legend('show');

% examine cdf of each cell type
nexttile
for j = 1:length(cellType)
   idx = Tc.identity == cellType(j) & Tc.selectivityPVal < pThr;

    % calculate cdf of effect sizes for each cell type, 
    [cdfEffectSize, xValuesCDF] = ecdf(Tc.selectivityEffectSize(idx));
    plot(xValuesCDF, cdfEffectSize, 'DisplayName', ['CDF of ' char(cellType(j))], 'LineWidth', 2);
    hold on;
end
xlabel('Selectivity Effect Size');
ylabel('Cumulative probability');
title('CDF of Effect Sizes by Cell Type');
legend('show');


%% perform population decoding                                                                                                                                                                                                        

nIter = 200; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
targetAmp = 5; % 0.33,0.76,1.6,3.5,7.7
usePCA = true; % use PCA step prior to model training
kStart = 1; % first PC to include
kComp = 15; % number of PCs to use in model training / testing
useFreqEffNeu = false; % use neurons with significant anova freq effects
useAllNeu = false; % use all available neurons for modeling (no random cell draws)
cellType = ["PV","SOM","EXC"]; % order for plotting

if multiAmp == true
    nCellDraws = 20:20:1000;
elseif multiAmp == false
    nCellDraws = 15:15:330;
end

% CREATE OUTPUT ARRAYS
% Performance of Cell Types (all)
% PCTall = []; %for populating with accuracy measurements
PCTall = zeros(nIter,nCellType,length(nCellDraws)); % (nj,ni,nd)
% Models (all)
MDLall = cell(nIter,nCellType,length(nCellDraws)); % (nj,ni,nd)

for nd = 1 : length(nCellDraws)
    disp(strcat("Draw ",num2str(nCellDraws(nd))))

    nCellDraw = nCellDraws(nd);
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
                idx = find(Tc.identity == cellType(ni) & Tc.selectivityPVal < pValThresh);
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

nCellDraws = 15:15:330;
figure('Color',[1 1 1],'Theme','Light'); hold on;
% plot fake data for correct legend order
plot([-1 -1],[-1 -1],'k-'); 
plot([-1 -1],[-1 -1],'k--');
fill([-1 -1 -1], [-1 -1 -1],colors.PV,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.SOM,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.EXC,'FaceAlpha',0.5,'EdgeColor','none');
for nt = 1 : size(PCTall,4)
    for nc = 1 : length(cellType)
        yMeans = mean(squeeze(PCTall(:,nc,:,nt)),1).*100;
        ySems = sem(squeeze(PCTall(:,nc,:,nt)),1).*100;
        xLooped = [nCellDraws flip(nCellDraws)]; 
        yLooped = [yMeans+ySems flip(yMeans-ySems)]; 
        fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.5,'EdgeColor','none');
        switch nt 
            case 1
                plot(nCellDraws,yMeans,'-','LineWidth',1,'Color',colors.(cellType(nc)));
            case 2
                plot(nCellDraws,yMeans,'--','LineWidth',1,'Color',colors.(cellType(nc)));
        end
    end
end
% plot([15 500],[20 20],'r--'); % chance performance level 
set(gca,'XLim',[15 330],'XTick',[15 75 125 175 225 275 330],...
    'YLim',[30 100],'YTick',30:10:100);
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Neurons used','FontName','Arial');
title('Frequency decoding performance by cell type','FontName','Arial');
subtitle('From neurons with frequency selectivity','FontName','Arial');
legend({"Freq. Mod.","All neurons","PV","SOM","EXC"},'Location','NorthWest');



























