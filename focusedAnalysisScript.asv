

%% load in datasets

% CHANGE DIRECTORIES (USER ADJUST)
cd('D:\Tactile_Synchrony\2P-Data\Cleaned\datasets');


% LOAD DATASET 1
% tic; load('dataset1_lightweight.mat'); toc;
% tic; load('dataset1_updated_20260612.mat'); toc;
tic; load('dataset1_updated_20260615.mat'); toc;

% LOAD DATASET 2
% tic; load("dataset2_lightweight.mat"); toc;
% tic; load("dataset2_lightweight_final.mat"); toc;
% tic; load("dataset2_updated_20260202.mat"); toc;
% tic; load("dataset2_updated_20260511.mat"); toc;
tic; load("dataset2_updated_20260601.mat"); toc;

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

% ESTABLISH AMPLITUDE VALUES FOR DATASET1
ampValsOg = [20; 35; 70; 45; 30];
d1AmpVals = zscore(ampValsOg);

% % PLOT
% figure('theme','light','color',[1 1 1]);
% plot(freqVal,ampValsOg,'ko-');
% set(gca,'Box','off','XLim',[250 1150],'YLim',[18 72],'XTick',freqVal,'YTick',[20:10:70],'FontName','Arial');
% ylabel('Amplitude (\mum)','FontName','Arial');
% xlabel('Frequency (Hz)','FontName','Arial');
% title('Stimulus Properties (dataset1)','FontName','Arial');


%% DECLARE ACTIVE DATASET

% USER INPUTS (CHANGE THESE)
Tc = T2; % <---- % T1 (Dataset1), T2 (dataset2)
minTrialNum = 15; % 25 (Dataset1) | 15 (Dataset2)

clear T1 T2 % remove others from memory to save RAM
% cellType = unique(Tc.identity);
cellType = ["PV","SOM","EXC"];
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
    freqVal = freqValAll(2:end); % skip 100HZ
    nFreq = length(freqVal);
end

if ~any(strcmp(fieldnames(Tc),'isTrained'))
    Tc.isTrained = Tc.trained;
    Tc.trained = [];
end

disp("Dataset Initiated")


% filter data 

minTrialCountsF = zeros(size(Tc,1),1); %for frequency alone
minTrialCounts = zeros(size(Tc,1),1); %for freq X amp

for i = 1:size(Tc,1)
    % Extract values
    cellFreq = Tc.frequencies{i};
    cellResp = Tc.responses{i};
    cellLateResp = Tc.lateResps{i};
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
            Tc.lateResps{i} = [];
            Tc.preStimVal{i} = [];
            Tc.amplitudes{i} = [];
        end

    else % DATASET1
        trialCountsF = accumarray(locF,1,[length(freqValAll),1]);
        % minimum trial number requirement for 300-1100 Hz
        minTrialCountsF(i) = min(trialCountsF(2:end));
        if  minTrialCountsF(i) >= minTrialNum && Tc.isTrained(i) == false
            % exclude 100 Hz
            Tc.frequencies{i} = cellFreq(locF>1); 
            Tc.responses{i} = cellResp(locF>1);
            Tc.lateResps{i} = cellLateResp(locF>1);
        else
            Tc.frequencies{i} = [];
            Tc.responses{i} = [];
            Tc.lateResps{i} = [];
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
is.trained = Tc.isTrained;

disp('cell identity');
tabulate(Tc.identity);



% calculate selectivity metric: anova effect size
% Pre-allocate columns in the table 

pValThresh = 0.05; % P VALUE THRESHOLD

Tc.responsivityEffectSize = zeros(nCellTotal, 1);
Tc.responsivityPVal = zeros(nCellTotal, 1);
corrValues = zeros(nCellTotal,1);

if multiAmp == true
    Tc.anovaFreqP = zeros(nCellTotal,1);
    Tc.anovaAmpP = zeros(nCellTotal,1);
    Tc.anovaMixP = zeros(nCellTotal,1);
    % ADD STRUCT FIELDS TO 'IS'
    is.FreqMod = false(nCellTotal,1); % frequency modulated ONLY
    is.AmpMod = false(nCellTotal,1); % amplitude modulated ONLY
    is.DualMod = false(nCellTotal,1); % Frequency AND Amplitude but NO interaction
    is.InterMod = false(nCellTotal,1); % Interaction effect
    is.nonMod = false(nCellTotal,1); % no effects
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
        if pVal(3) < pValThresh
            is.InterMod(i) = true;
        elseif pVal(1) < pValThresh && pVal(2) > pValThresh
            is.FreqMod(i) = true;
        elseif pVal(1) > pValThresh && pVal(2) < pValThresh
            is.AmpMod(i) = true;
        elseif pVal(1) < pValThresh && pVal(2) < pValThresh
            is.DualMod(i) = true;
        else
            is.nonMod(i) = true;
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

        % correlate with mean response per freq known amp values per freq
        meanResp = zeros(nFreq,1);
        for nf = 1 : nFreq
            meanResp(nf) = mean(cellResp(cellFreq == freqVal(nf)));
        end
        % corrValues(i) = corr(d1AmpVals,zscore(meanResp)); %z-score values then pearson corr
        corrValues(i) = corr(meanResp,ampValsOg,"Type","Spearman"); 
    
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

% calculate custom selectivity metric for each neuron in Tc
avgRespAll = zeros(size(Tc,1),length(freqVal));
prefFreqAll = zeros(size(Tc,1),1);
for n = 1 : size(Tc,1)
    locResps = Tc.responses{n};
    locFreqs = Tc.frequencies{n};
    avgRespAll(n,:) = zeros(length(freqVal),1); % average dF/F per frequency
    for nf = 1 : length(freqVal)
        avgRespAll(n,nf) = mean(locResps(locFreqs == freqVal(nf)));
    end
    % Tc.selectivityMetric(n) = quickSelectivity(avgRespAll(n,:)); % calculate selectivity metric
    Tc.selectivityMetric(n) = quickSelectivityEntr(avgRespAll(n,:)); % calculate selectivity metric
    [~,prefFreqAll(n)] = max(avgRespAll(n,:));
end

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
disp('Selectivity & ANOVAs calculated'); fprintf('\n');


if multiAmp == 0
for nm = 1 : size(is.mouseNum,2)
    disp(mouseNums(nm));
    nPV = sum(is.mouseNum(:,nm) & is.PV);
    nSelPV = sum(is.mouseNum(:,nm) & is.PV & is.selective);
    disp(strcat(num2str(nPV)," PV, ",num2str(nSelPV)," sel. PV"));
    nSOM = sum(is.mouseNum(:,nm) & is.SOM);
    nSelSOM = sum(is.mouseNum(:,nm) & is.SOM & is.selective);
    disp(strcat(num2str(nSOM)," SOM, ",num2str(nSelSOM)," sel. SOM"));
    nEXC = sum(is.mouseNum(:,nm) & is.EXC);
    nSelEXC = sum(is.mouseNum(:,nm) & is.EXC & is.selective);
    disp(strcat(num2str(nEXC)," EXC, ",num2str(nSelEXC)," sel. EXC"));
    if nm ~= size(is.mouseNum,2)
        disp('------------------------');
    end
end
end


 %% LOAD TO HERE

% cd('C:\Users\skich\Desktop\WORK');
% cd('C:\Users\skich\Box\Tactile_Synchrony\2P-Data\Cleaned\datasets');
cd('D:\Tactile_Synchrony\2P-Data\Cleaned\datasets')

load('dataset2_forFigures.mat');

% load('dataset1_forFigures.mat');


% addpath('C:\Users\skich\Box\Tactile_Synchrony\Haptics-2P\Analysis');
% cd('C:\Users\skich\Documents\GitHub\neural-analysis-suite');

% % ANALYSIS SCRIPS
%
% FAS_data1_removeAmpCorr.m   
% FAS_multiCompare.m    
% FAS_simultaneousDecode.m   
% FAS_decodeAmp.m   
% FAS_omit700Hz.m     
% FAS_titrateNumPCs.m
% FAS_decodeFreqWithAmp.m  
% FAS_selTypePerf_allCells.m  
% FAS_titrateSelectivity.m 
% FAS_indvMouseOmit.m  
% FAS_selTypePerformance.m 
% FAS_tuningHeatmap.m 
% FAS_Dataset2Figures.m
% FAS_trainingEffectFigures.m
% FAS_selTypePerf_allCells_ampAndFreq.m


%% NEURON AMPLITUDE CORRELATIONS PLOT

figure('theme','light','color',[1 1 1]); tiledlayout(3,1);
nexttile
h1 = histogram(corrValues(is.PV & is.selective),'Normalization','probability','FaceColor',colors.PV,'BinWidth',0.02);
set(gca,"YLim",[0 0.1],"YTick",0:0.02:0.1,"YTickLabel",0:2:10,'XLim',[0 1],'box','off','XTickLabel',[]);
title("Neurons with responses correlated to amplitude level",'FontName','Arial');
nexttile
h2 = histogram(corrValues(is.SOM & is.selective),'Normalization','probability','FaceColor',colors.SOM,'BinWidth',0.02);
set(gca,"YLim",[0 0.1],"YTick",0:0.02:0.1,"YTickLabel",0:2:10,'XLim',[0 1],'box','off','XTickLabel',[]);
ylabel("Percentage of neurons (%)",'FontName','Arial');
nexttile
h3 = histogram(corrValues(is.EXC & is.selective),'Normalization','probability','FaceColor',colors.EXC,'BinWidth',0.02);
set(gca,"YLim",[0 0.1],"YTick",0:0.02:0.1,"YTickLabel",0:2:10,'XLim',[0 1],'box','off');
xlabel("Correlation value (r-squared)",'FontName','Arial');

figure('theme','light','color',[1 1 1]);
tiledlayout(1,3);
nexttile;
plot(corrValues(is.PV & is.selective),Tc.selectivityEffectSize(is.PV & is.selective),'.','Color',colors.PV);
set(gca,'XLim',[0 1],'YLim',[0 1]);
nexttile;
plot(corrValues(is.SOM & is.selective),Tc.selectivityEffectSize(is.SOM & is.selective),'.','Color',colors.SOM);
set(gca,'XLim',[0 1],'YLim',[0 1]);
nexttile;
plot(corrValues(is.EXC & is.selective),Tc.selectivityEffectSize(is.EXC & is.selective),'.','Color',colors.EXC);
set(gca,'XLim',[0 1],'YLim',[0 1]);

%% BAR PLOTS: RESPONSIVE(%) and SELECTIVE(%) FROM CELL TYPES

if multiAmp == true
    Tc.selectivityPVal = Tc.anovanFreqP;
end
% % PLOT TOGETHER IN ONE FIGURE
% figure('theme','light','color',[1 1 1]); hold on;
% respPerc = zeros(length(cellType),1);
% selPerc = zeros(length(cellType),1);
% xVals = [1; 2; 3; 5; 6; 7];
% for i = 1:length(cellType)
%     % RESPONSIVE UNITS CALC
%     indx = Tc.identity == cellType(i); % target neurons
%     numResps = sum(Tc.responsivityPVal(indx)<0.05); % number of responsive neurons
%     respPerc(i) = numResps/sum(indx)*100; % responsive neuron percentage
%     respSTD = mannyPropSTD(numResps,sum(indx)) * 100; % responsive neuron std
%     % RESPOINSIVE UNITS PLOT
%     bar(xVals(i),respPerc(i),'FaceColor',colors.(cellType(i))); 
%     text(xVals(i),(respPerc(i)+respSTD+0.02),strcat("n=",num2str(numResps)),'HorizontalAlignment',...
%         'center','VerticalAlignment','bottom','FontName','Arial','FontSize',9)
%     errorbar(xVals(i),respPerc(i),respSTD,Color=[0,0,0],LineWidth=1.4);
%     % SELECTIVE UNITS CALC
%     numSels = sum(Tc.selectivityPVal(indx)<0.05); % number of selective neurons
%     selPerc(i) = numSels/numResps*100;
%     selSTD = mannyPropSTD(numSels,numResps) * 100; % selective neuron std
%     % SELECTIVE UNITS PLOT
%     bar(xVals(i+3),selPerc(i),'FaceColor',colors.(cellType(i))); 
%     text(xVals(i+3),(selPerc(i)+selSTD+0.02),strcat("n=",num2str(numSels)),'HorizontalAlignment',...
%         'center','VerticalAlignment','bottom','FontName','Arial','FontSize',9)
%     errorbar(xVals(i+3),selPerc(i),selSTD,Color=[0,0,0],LineWidth=1.4);
% end
% set(gca,'Xtick',[2 6],'XTickLabels',{"Responsive","Selective"},'ylim',[20 100],'FontName','Arial','Box','off');
% title('Response properties of all neurons','FontName','Arial');
% ylabel('percentage of neurons per group','FontName','Arial');
% legend({'PV','','','','SOM','','','','EXC'},'Location','NorthEast','box','off');

% PLOT RESPONSIVITY
figure('theme','light','color',[1 1 1]); hold on;
respPerc = zeros(length(cellType),1);
for i = 1:length(cellType)
    % RESPONSIVE UNITS CALC
    indx = Tc.identity == cellType(i); % target neurons
    numResps = sum(Tc.responsivityPVal(indx)<0.05); % number of responsive neurons
    respPerc(i) = numResps/sum(indx); % responsive neuron percentage
    respSE = mannyPropSTD(numResps,sum(indx)); % responsive neuron estimated standard error
    respCI = respSE * 1.96; %SE to 95% CI
    % RESPOINSIVE UNITS PLOT
    bar(i,respPerc(i)*100,'FaceColor',colors.(cellType(i))); 
    text(i,((respPerc(i)*100)+(respCI*100)+0.02),strcat("n=",num2str(numResps)),'HorizontalAlignment',...
        'center','VerticalAlignment','bottom','FontName','Arial','FontSize',9)
    errorbar(i,respPerc(i)*100,respCI*100,Color=[0,0,0],LineWidth=1.4);
end
set(gca,'XTick',[1 2 3],'XTickLabels',cellType,'ylim',[70 100],'FontName','Arial','Box','off');
title('Neurons with stimulus response','FontName','Arial');
ylabel('Percentage of neurons per group','FontName','Arial');


% PLOT SELECTIVITY
figure('theme','light','color',[1 1 1]); hold on;
selPerc = zeros(length(cellType),1);
for i = 1:length(cellType)
    % SELECTIVE UNITS CALC
    indx = Tc.identity == cellType(i); % target neurons
    numResps = sum(Tc.responsivityPVal(indx)<0.05); % number of responsive neurons
    numSels = sum(Tc.selectivityPVal(indx)<0.05); % number of selective neurons
    selPerc(i) = numSels/numResps;
    selSE = mannyPropSTD(numSels,numResps); % selective neuron estimated standard error
    selCI = selSE * 1.96;%SE to 95% CI
    % SELECTIVE UNITS PLOT
    bar(i,selPerc(i)*100,'FaceColor',colors.(cellType(i))); 
    text(i,((selPerc(i)*100)+(selCI*100)+0.02),strcat("n=",num2str(numSels)),'HorizontalAlignment',...
        'center','VerticalAlignment','bottom','FontName','Arial','FontSize',9)
    errorbar(i,selPerc(i)*100,selCI*100,Color=[0,0,0],LineWidth=1.4);
end
set(gca,'Xtick',[1 2 3],'XTickLabels',cellType,'ylim',[20 60],'FontName','Arial','Box','off');
title('Neurons with frequency selectivity','FontName','Arial');
ylabel('Percentage of responsive neurons','FontName','Arial');


%% VENN DIAGRAM PER CELL

nc = 3;

nF = sum(is.(cellType(nc)) & Tc.anovanFreqP<0.05);
nA = sum(is.(cellType(nc)) & Tc.anovanAmpP<0.05);
nI = sum(is.(cellType(nc)) & Tc.anovanMixP<0.05);

nFA = sum(is.(cellType(nc)) & Tc.anovanFreqP<0.05 & Tc.anovanAmpP<0.05 & Tc.anovanMixP>0.05);
nFI = sum(is.(cellType(nc)) & Tc.anovanFreqP<0.05 & Tc.anovanAmpP>0.05 & Tc.anovanMixP<0.05);
nAI = sum(is.(cellType(nc)) & Tc.anovanFreqP>0.05 & Tc.anovanAmpP<0.05 & Tc.anovanMixP<0.05);
nFAI = sum(is.(cellType(nc)) & Tc.anovanFreqP<0.05 & Tc.anovanAmpP<0.05 & Tc.anovanMixP<0.05);

nFo = nF - nFA - nFI - nFAI;
nAo = nA - nFA - nAI - nFAI;
nIo = nI - nFI - nAI - nFAI;

% A = [nF nA nI] % F, A, F*A
% I = [nFA nFI nAI nFAI] % interactions
A = [40 40 30];
I = [20 5 5 15];

figure('theme','light','color',[1 1 1],'Position',[100 100 550 500]);
venn(A,I,'FaceAlpha',0.5,'FaceColor',colors.(cellType(nc)));
set(gca,'XColor',[1 1 1],'YColor',[1 1 1]);

annotation("textbox", [0.4503 0.343 0.07063 0.05], "String", num2str(nFA), "FontName", "Arial", "HorizontalAlignment", "center", "LineStyle", "none")
annotation("textbox", [0.2593 0.203 0.06599 0.054], "String", "F", "FontName", "Arial", "FontSize", 12, "FontWeight", "bold", "HorizontalAlignment", "center", "LineStyle", "none")
annotation("textbox", [0.2356 0.381 0.05762 0.05], "String", num2str(nFo), "FontName", "Arial", "HorizontalAlignment", "center", "LineStyle", "none")
annotation("textbox", [0.6761 0.361 0.05762 0.05], "String", num2str(nAo), "FontName", "Arial", "HorizontalAlignment", "center", "LineStyle", "none")
annotation("textbox", [0.4359 0.641 0.09944 0.054], "String", "F*A", "FontName", "Arial", "FontSize", 12, "FontWeight", "bold", "HorizontalAlignment", "center", "LineStyle", "none")
annotation("textbox", [0.6506 0.199 0.07156 0.054], "String", "A", "FontName", "Arial", "FontSize", 12, "FontWeight", "bold", "HorizontalAlignment", "center", "LineStyle", "none")
annotation("textbox", [0.4303 0.211 0.1069 0.054], "String", "F+A", "FontName", "Arial", "FontSize", 12, "FontWeight", "bold", "HorizontalAlignment", "center", "LineStyle", "none")
annotation("textbox", [0.441 0.517 0.0855 0.05], "String", num2str(nFAI), "FontName", "Arial", "HorizontalAlignment", "center", "LineStyle", "none")
annotation("textbox", [0.4531 0.723 0.05762 0.05], "String", num2str(nIo), "FontName", "Arial", "HorizontalAlignment", "center", "LineStyle", "none")
annotation("textbox", [0.5488 0.555 0.07063 0.05], "String", num2str(nAI), "FontName", "Arial", "HorizontalAlignment", "center", "LineStyle", "none")
annotation("textbox", [0.3388 0.555 0.07063 0.05], "String", num2str(nFI), "FontName", "Arial", "HorizontalAlignment", "center", "LineStyle", "none")
hTextboxshape = findall(gcf,"Type","textboxshape");






%% SELECTIVITY METRIC PROBABILITY DENSITY FUNCTION
% examine pdf of each cell type
figure('theme','light','color',[1 1 1]); hold on;
pThr = 0.05;
for nc = 1:length(cellType)
    idx = Tc.identity == cellType(nc) & is.selective;
    [pdfEffectSize, xValues] = ksdensity(Tc.selectivityMetric(idx),'Support',[0 1]);
    plot(xValues, pdfEffectSize,'LineWidth',2,'Color',colors.(cellType(nc)));
end
xlabel('Selectivity Strength');
ylabel('Probability Density');
title('PDF of Effect Sizes by Cell Type','FontName','Arial');
legend({'PV','SOM','EXC'},'Location','northwest','box','off'); 
hold off;


%% SELECTIVITY METRIC CUMMULATIVE DISTRIBUTION FUNCTION
figure('theme','light','color',[1 1 1]); hold on;
plot([0 1],[0.33 0.33],'LineStyle',':','Color',repmat(0.7,1,3),'LineWidth',0.5);
plot([0 1],[0.66 0.66],'LineStyle',':','Color',repmat(0.7,1,3),'LineWidth',0.5);
text(0.8,0.05,'bottom third','HorizontalAlignment','right','FontName','Arial','Color',repmat(0.7,1,3));
text(0.8,0.38,'middle third','HorizontalAlignment','right','FontName','Arial','Color',repmat(0.7,1,3));
text(0.8,0.71,'top third','HorizontalAlignment','right','FontName','Arial','Color',repmat(0.7,1,3));
for nc = 1:length(cellType)
    idx = Tc.identity == cellType(nc) & is.selective;
    [cdfEffectSize,xValuesCDF] = ecdf(Tc.selectivityMetric(idx));
    plot(xValuesCDF,cdfEffectSize,'LineWidth',2,'Color',colors.(cellType(nc)));
end
xlabel('Selectivity Strength','FontName','Arial');
ylabel('Cumulative probability','FontName','Arial');
title(["Frequency selectivity strength distributions", "by cell type"],'FontName','Arial');
legend({'','','PV','SOM','EXC'},'Location','northeast','box','off','FontName','Arial'); 
set(gca,'YLim',[0 1.01],'XLim',[0.1 0.8]);
hold off;

%% PERMUTATION TEST BETWEEN SEL POPS
allSelPV = Tc.selectivityMetric(Tc.identity == "PV" & is.selective);
allSelSOM = Tc.selectivityMetric(Tc.identity == "SOM" & is.selective);
allSelEXC = Tc.selectivityMetric(Tc.identity == "EXC" & is.selective);
figure('Color',[1 1 1]); hold on;
swarmchart(ones(size(allSelPV)),allSelPV,'filled','o','MarkerFaceColor',colors.PV,...
    'MarkerEdgeAlpha',0,'SizeData',4);
swarmchart(ones(size(allSelSOM))*2,allSelSOM,'filled','o','MarkerFaceColor',colors.SOM,...
    'MarkerEdgeAlpha',0,'SizeData',4);
swarmchart(ones(size(allSelEXC))*3,allSelEXC,'filled','o','MarkerFaceColor',colors.EXC,...
    'MarkerEdgeAlpha',0,'SizeData',4,'MarkerFaceAlpha',0.3);
plot([1.1 1.9],[0.55 0.55],'k-','LineWidth',1); 
text(1.5,0.55,'n.s.','HorizontalAlignment','center','VerticalAlignment','bottom');
plot([2.1 2.9],[0.55 0.55],'k-','LineWidth',1);
text(2.5,0.55,'n.s.','HorizontalAlignment','center','VerticalAlignment','bottom');
plot([1.1 2.9],[0.6 0.6],'k-','LineWidth',1);
text(2,0.6,'n.s.','HorizontalAlignment','center','VerticalAlignment','bottom');
set(gca,'XTick',1:3,'XTickLabel',cellType,'XLim',[0.5 3.5]);
ylabel('Selectivity Strength','FontName','Arial');
title(["Frequency selectivity strength", "by cell type"],'FontName','Arial');

%% EXAMPLE SEL VALUED NEURONS

figure('Color',[1 1 1]); hold on;
xi = [1 2 3; 4 5 6];
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
targ = [0.4 0.15];
for nc = 1 : length(cellType)
    for nt = 1 : length(targ)
        nexttile(xi(nt,nc)); hold on;
        selNs = find(Tc.identity == cellType(nc) & is.selective);
        [val1,idx1] = sort(Tc.selectivityMetric(selNs),'descend');
        [~,idx2] = min(abs(targ(nt) - val1));
        exResps = Tc.responses{selNs(idx1(idx2))}; % close to 0.4
        exFreqs = Tc.frequencies{selNs(idx1(idx2))}; % close to 0.4
        for nf = 1 : length(freqVal)
            locResps = exResps(exFreqs == freqVal(nf));
            locMean = mean(locResps);
            locCI = locMean - mkCI(locResps);
            plot(freqVal(nf),locMean,'.','MarkerSize',8,'Color',colors.(cellType(nc)))
            errorbar(freqVal(nf),locMean,locCI(1),'vertical','Color',colors.(cellType(nc)));
        end
        if nt == 2
            set(gca,'XLim',[200 1200],'YLim',[-0.02 0.14],'XTick',freqVal);
        else
            set(gca,'XLim',[200 1200],'YLim',[-0.02 0.14],'XTick',freqVal,'XTickLabels',[]);
        end
        if nc == 1
            ylabel("Frequency (Hz)",'FontName','Arial');
        elseif nc == 2
            title(strcat("Selectivity strength = ",num2str(targ(nt))),'FontName','Arial');
            if nt == 2
                xlabel("Frequency (Hz)",'FontName','Arial');
            end
        end
    end
end


%% INSET PLOTS OF POPULATION THIRDS

nBins = 3; % number of bins to divide into
meanResps = struct();
meanCIs = struct();
bins = cell(nBins,length(cellType));
figure('theme','light','color',[1 1 1]);
tiledlayout(nBins,1);
plotOrder = flip(1:nBins); % plot reverse order: top third, middle third, bottom third

for nc = 1:length(cellType)
    idx = Tc.identity == cellType(nc) & is.selective;
    sortedSelVals = sort(Tc.selectivityMetric(idx));
    minBinSize = floor(length(sortedSelVals)/nBins); % rough bin size
    mod1 = mod(length(sortedSelVals),nBins); % remainder
    binSizes = [repmat(minBinSize,1,nBins)]; % distribute bins
    binSizes(1:mod1) = binSizes(1:mod1) + 1; % add in remainder(s)
    bins(:,nc) = mat2cell(sortedSelVals, binSizes, 1);
    % calculate means with confidence intervals
    meanResps.(cellType(nc)) = cellfun(@mean,bins(:,nc));
    meanCIs.(cellType(nc)) = reshape(cell2mat(cellfun(@mkCI,bins(:,nc),'UniformOutput',false)),...
        2,nBins)';
end

for ni = 1 : nBins
    lni = plotOrder(ni);
    yVals = [];
    nexttile(ni); hold on;
    for nc = 1 : length(cellType)
        locErr = meanCIs.(cellType(nc))(lni,:) - meanResps.(cellType(nc))(lni);
        errorbar(nc,meanResps.(cellType(nc))(lni),locErr(1),locErr(2),'vertical',...
            'Color',colors.(cellType(nc)),'LineWidth',1);
        plot(nc,meanResps.(cellType(nc))(lni),'.','Color',colors.(cellType(nc)),...
            'MarkerSize',10);
        yVals = [yVals meanResps.(cellType(nc))(lni)+locErr];
    end
    % fix labels
    if ni == 1
        title('Top third','FontName','Arial');
        yTL = string(num2str([min(yVals); max(yVals)],2));
        set(gca,'XLim',[0.5 3.5],'YTick',[min(yVals) max(yVals)],'YTickMode','Manual',...
            'XColor','none','YLim',[min(yVals)-.005 max(yVals)+.005],'YTickLabel',yTL,...
            'YTickLabelRotation',0,'FontSize',10);
    elseif ni == 2
        title('Middle third','FontName','Arial');
        ylabel('Selectivity strength','FontName','Arial');
        yTL = string(num2str([min(yVals); max(yVals)],2));
        set(gca,'XLim',[0.5 3.5],'YTick',[min(yVals) max(yVals)],'YTickMode','Manual',...
            'XColor','none','YLim',[min(yVals)-.005 max(yVals)+.005],'YTickLabel',yTL,...
            'YTickLabelRotation',0,'FontSize',10);
    elseif ni == 3
        title('Bottom third','FontName','Arial');
        yTL = string(num2str([min(yVals); max(yVals)],2));
        set(gca,'XLim',[0.5 3.5],'YTick',[min(yVals) max(yVals)],'YTickMode','Manual',...
            'XTick',[1 2 3],'XTickLabels',cellType,'YLim',[min(yVals)-.005 max(yVals)+.005],...
            'YTickLabel',yTL,'YTickLabelRotation',0,'FontSize',10);
    end

    % stats PV vs SOM
    [p1,h1] = ranksum(bins{lni,1},bins{lni,2});
    if h1 == 1
        % maxY = max([meanCIs.PV(lni,:) meanCIs.SOM(lni,:)]) + 0.001;
        maxY = max([meanCIs.PV(lni,:) meanCIs.SOM(lni,:)]);
        plot([1.1 1.9],[maxY maxY],'k-');
        text(1.5,maxY,'*','HorizontalAlignment','center','FontSize',12);
    end
    % stats PV vs EXC
    [p2,h2] = ranksum(bins{lni,1},bins{lni,3});
    if h2 == 1
        % maxY = max([meanCIs.PV(lni,:) meanCIs.EXC(lni,:)]) + (0.0015*lni);
        maxY = max([meanCIs.PV(lni,:) meanCIs.EXC(lni,:)])
        plot([1.1 2.9],[maxY maxY],'k-');
        text(2,maxY,'*','HorizontalAlignment','center','FontSize',12);
    end
    % stats SOM vs EXC
    [p3,h3] = ranksum(bins{lni,2},bins{lni,3});
    if h3 == 1
        % maxY = max([meanCIs.SOM(lni,:) meanCIs.EXC(lni,:)]) + 0.001;
        maxY = max([meanCIs.SOM(lni,:) meanCIs.EXC(lni,:)]) - 0.002*lni;
        plot([2.1 2.9],[maxY maxY],'k-');
        text(2.5,maxY,'*','HorizontalAlignment','center','FontSize',12);
    end
end


%% SELECTIVITY METRIC AS PDF BY EVEN BINS

nBins = 3; % number of bins to divide into
pThr = 0.05; % p valuse threshold
nBoot = 1000; % number of bootstrap iterations
xGrid = linspace(0,1,100); % fixed evaluation points for ksdensity
nTypes = numel(cellType);
figure('theme','light','color',[1 1 1]); hold on;
tiledlayout(3,1);
for nc = 1:length(cellType)
    idx = Tc.identity == cellType(nc) & is.selective;
    sortedSelVals = sort(Tc.selectivityMetric(idx));
    minBinSize = floor(length(sortedSelVals)/nBins); % rough bin size
    mod1 = mod(length(sortedSelVals),nBins); % remainder
    binSizes = [repmat(minBinSize,1,nBins)]; % distribute bins
    binSizes(1:mod1) = binSizes(1:mod1) + 1; % add in remainder(s)
    bins = mat2cell(sortedSelVals, binSizes, 1); % populate bins
    for nb = 1 : nBins
        locVals = bins{nb};
        % Compute observed PDF 
        [pdfMean, ~] = ksdensity(locVals, xGrid, 'Support', [0 1]);
        % Bootstrap resampling
        bootPDFs = zeros(nBoot, numel(xGrid));
        switch nc
            case 1
                bootPDFs_PV = bootPDFs;
            case 2
                bootPDFs_SOM = bootPDFs;
            case 3
                bootPDFs_EXC = bootPDFs;
        end
        nVals = numel(locVals);
        for b = 1:nBoot
            resampleVals = locVals(randi(nVals, [nVals, 1])); % sample with replacement
            bootPDFs(b,:) = ksdensity(resampleVals, xGrid, 'Support', [0 1]);
        end
        % Compute confidence intervals
        pdfLow  = prctile(bootPDFs, 2.5, 1);
        pdfHigh = prctile(bootPDFs, 97.5, 1);
        nexttile(nb); hold on;
        % Plot mean curve
        plot(xGrid, pdfMean, 'LineWidth', 2, 'Color', colors.(cellType(nc)), ...
             'DisplayName', char(cellType(nc)));
        % Plot shaded CI region
        fill([xGrid fliplr(xGrid)],[pdfLow fliplr(pdfHigh)],colors.(cellType(nc)), ...
            'FaceAlpha',0.2,'EdgeColor','none','DisplayName','95% CI');
        hold off;
    end
end

% set(gca,'XTick',xMat(:,2),'YLim',[0.12 0.52],'XTickLabel',{'Q1','Q2','Q3','Q4','Q5'});
% xlabel('Quintile bin','FontName','Arial');
% ylabel('selectivity strength','FontName','Arial');
% title('Frequency Selectivity by quartile','FontName','Arial');
% % legend('show','Location','East','box','off','FontName','Arial'); 
% hold off;


%% SELECTIVITY METRIC SPLIT BY EVEN BINS

xSp = 10; % spacing for plotting
nBins = 3; % number of bins to divide into

nTypes = numel(cellType);
d1 = nTypes + xSp;
xMat = (1:d1*nBins);
xMat = reshape(xMat, d1, nBins).';
xMat = xMat(:, 1:nTypes);
figure('theme','light','color',[1 1 1]); hold on;
for nc = 1:length(cellType)
    idx = Tc.identity == cellType(nc) & is.selective;
    sortedSelVals = sort(Tc.selectivityMetric(idx));
    minBinSize = floor(length(sortedSelVals)/nBins); % rough bin size
    mod1 = mod(length(sortedSelVals),nBins); % remainder
    binSizes = [repmat(minBinSize,1,nBins)]; % distribute bins
    binSizes(1:mod1) = binSizes(1:mod1) + 1; % add in remainder(s)
    bins = mat2cell(sortedSelVals, binSizes, 1); % populate bins
    % calculate means with confidence intervals
    meanResps = cellfun(@mean,bins);
    meanCIs = reshape(cell2mat(cellfun(@mkCI,bins,'UniformOutput',false)),2,nBins)';
    for ni = 1 : length(meanResps)
        plot( repmat(xMat(ni,nc),2,1),meanCIs(ni,:),'LineStyle','-','Color',...
            colors.(cellType(nc)),'LineWidth',1.2); % CI line
        plot(xMat(ni,nc),meanResps(ni),'.','Color',colors.(cellType(nc)),...
            'MarkerSize',8);
    end
    %plot(xMat(:,nc),meanResps,'-','Color',colors.(cellType(nc)),'LineWidth',0.5);
end
set(gca,'XTick',xMat(:,2),'YLim',[0.12 0.52],'XTickLabel',{'Q1','Q2','Q3','Q4','Q5'});
xlabel('Quintile bin','FontName','Arial');
ylabel('selectivity strength','FontName','Arial');
title('Frequency Selectivity by quartile','FontName','Arial');
% legend('show','Location','East','box','off','FontName','Arial'); 
hold off;


%% PDF WITH PERMUTATIONS 
pThr = 0.05; % p valuse threshold
nBoot = 1000; % number of bootstrap iterations
xGrid = linspace(0,1,100); % fixed evaluation points for ksdensity
figure('theme','light','color',[1 1 1]); hold on;

for nc = 1:length(cellType)
    % Extract data for this cell type
    idx = Tc.identity == cellType(nc) & Tc.selectivityPVal < pThr;
    locVals = Tc.selectivityMetric(idx);
    % Compute observed PDF 
    [pdfMean, ~] = ksdensity(locVals, xGrid, 'Support', [0 1]);
    % Bootstrap resampling
    bootPDFs = zeros(nBoot, numel(xGrid));
    switch nc
        case 1
            bootPDFs_PV = bootPDFs;
        case 2
            bootPDFs_SOM = bootPDFs;
        case 3
            bootPDFs_EXC = bootPDFs;
    end
    nVals = numel(locVals);
    for b = 1:nBoot
        resampleVals = locVals(randi(nVals, [nVals, 1])); % sample with replacement
        bootPDFs(b,:) = ksdensity(resampleVals, xGrid, 'Support', [0 1]);
    end
    % Compute confidence intervals
    pdfLow  = prctile(bootPDFs, 2.5, 1);
    pdfHigh = prctile(bootPDFs, 97.5, 1);
    % Plot mean curve
    plot(xGrid, pdfMean, 'LineWidth', 2, 'Color', colors.(cellType(nc)), ...
         'DisplayName', char(cellType(nc)));
    % Plot shaded CI region
    fill([xGrid fliplr(xGrid)],[pdfLow fliplr(pdfHigh)],colors.(cellType(nc)), ...
        'FaceAlpha',0.2,'EdgeColor','none','DisplayName','95% CI');
end

xlabel('Selectivity Strength','FontName','Arial');
ylabel('Probability Density','FontName','Arial');
legend({'PV','','SOM','','EXC','95% CI'},'Location','East','box','off','FontName','Arial'); 
title('Frequency selectivity strength of each cell type','FontName','Arial');


%% PDF PERMUTATION INSET PLOT FOR SIGNIFICANCE
figure('theme','light','color',[1 1 1]);
tiledlayout(2,2);
uTile = [1 2 4];
for nc = 1 : 3
    if nc == 1 % PV vs SOM
        diffBoot = bootPDFs_PV - bootPDFs_SOM;  % difference in bootstrapped PDFs
    elseif nc == 2 % PV vs EXC
        diffBoot = bootPDFs_PV - bootPDFs_EXC;
    elseif nc == 3 % SOM vs EXC
        diffBoot = bootPDFs_SOM - bootPDFs_EXC;
    end
    diffLow  = prctile(diffBoot, 2.5, 1);
    diffHigh = prctile(diffBoot, 97.5, 1);
    diffMean = mean(diffBoot, 1);
    % "Significant divergence" where 95% CI does NOT cross zero
    sigIdx = diffLow > 0 | diffHigh < 0;
    nexttile(uTile(nc)); hold on;
    plot(xGrid, diffMean, 'k', 'LineWidth', 2);
    fill([xGrid fliplr(xGrid)], [diffLow fliplr(diffHigh)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    yline(0, '--');
    plot(xGrid(sigIdx), diffMean(sigIdx), 'r.', 'MarkerSize', 10);
    xlabel('Selectivity Metric','FontName','Arial');
    ylabel('PDF difference (SOM - EXC)','FontName','Arial');
    if nc == 1
        title({'Pairwise PDF divergence','PV vs. SOM'},'FontName','Arial');
    elseif nc == 2
        title('PV vs. EXC','FontName','Arial');
    elseif nc == 3
        title('SOM vs. EXC','FontName','Arial');
    end
end
nexttile(3); hold on;
plot([-1 -1],[0 0],'k-','LineWidth',2);
fill([-1 -1 -1 -1],[-1 -1 -1 -1],'k','FaceAlpha',0.2,'EdgeColor','none');
plot(-1,-1,'r.','MarkerSize', 10);
legend({'PDF diff.','95% CI','Significance'},'Location','NorthEast','box','off','FontName','Arial');
set(gca,'XLim',[1 2],'XColor','none','YColor','none');


%% Frequency tuning heatmaps

% calculate custom selectivity metric for each neuron in Tc
avgRespAll = zeros(size(Tc,1),length(freqVal));
prefFreqAll = zeros(size(Tc,1),1);
for n = 1 : size(Tc,1)
    locResps = Tc.responses{n};
    locFreqs = Tc.frequencies{n};
    avgRespAll(n,:) = zeros(length(freqVal),1); % average dF/F per frequency
    for nf = 1 : length(freqVal)
        avgRespAll(n,nf) = mean(locResps(locFreqs == freqVal(nf)));
    end
    Tc.selectivityMetric(n) = quickSelectivity(avgRespAll(n,:)); % calculate selectivity metric
    [~,prefFreqAll(n)] = max(avgRespAll(n,:));
end

figure('theme','light','color',[1 1 1]);
tiledlayout(1,3,'TileSpacing','tight','padding','tight');
for nc = 1 : length(cellType)
    locHeatmap = [];
    for nf = 1 : length(freqVal)
        indx = is.(cellType(nc)) & is.selective & prefFreqAll == nf;
        if sum(indx) == 0
            continue
        end
        locMat = avgRespAll(indx,:);
        [~,sIndx] = sort(Tc.selectivityMetric(indx));
        sortedMat = normalize(locMat(sIndx,:),2,'range');
        locHeatmap = cat(1,locHeatmap,sortedMat);
    end
    nexttile(nc); hold on;
    imagesc(flip(locHeatmap,1));
    subtitle(cellType(nc),'FontName','Arial');
    set(gca,'XTick',1:5,'XTickLabel',freqVal,'XTickLabelRotation',45,'YLim',...
        [0 size(locHeatmap,1)],'FontName','Arial');
    if nc == 1
        newYTicks = round(sum(is.(cellType(nc)) & is.selective) .* (0:0.2:1));
        set(gca,'YTick',newYTicks,'YTickLabel',flip(0:20:100),'FontName','Arial','FontSize',12);
        ylabel('Percentage of neurons (%)','FontName','Arial','FontSize',12);
    elseif nc == 2
        title('Sorted frequency tuning per cell type','FontName','Arial','FontSize',12);
        set(gca,'YColor','none','FontName','Arial','Fontsize',12);
        xlabel('Frequency (Hz)','FontName','Arial','FontSize',12)
    elseif nc == 3
        set(gca,'YColor','none','FontName','Arial','FontSize',12);
    end
end



%% FREQ PREFERENCE CHANGING WITH AMPLITUDE
nc = 3;
targets = find(is.DualMod & is.(cellType(nc))); % identify target cells
allSlopes = NaN(size(targets));
for nn = 1 : length(targets) 
    prefFreq = Tc.perfFreqPerAmp{targets(nn)}; % preferred freq at each amp
    if sum(prefFreq ~=0) >= 2 % check for at least 2 preferred freq values across amplitude
        trueFreq = NaN(5,1); % yVals (freq)
        for nf = 1 : length(prefFreq)
            if prefFreq(nf) == 0
                continue; %skip this amp
            end
            trueFreq(nf) = freqVal(prefFreq(nf));
        end
        tbl1 = table(ampValAll,trueFreq,'VariableNames',{'x','y'}); % fit mean vals (visualization only)
        f1 = fitlm(tbl1,'y ~ x'); % linear fit
        allSlopes(nn) = f1.Coefficients.Estimate(2); % slope
        % [yPred,yCI] = predict(f1,ampValAll,'Alpha',0.05);
        % plot(ampValAll,yPred,'--','LineWidth',0.5,'Color',colors.(cellType(nc))); % fitted line
    end
end

figure('theme','light','color',[1 1 1]); hold on; 
% histogram(allSlopes,'NumBins',100);
swarmchart(ones(size(allSlopes)),allSlopes)

targets = find(is.InterMod & is.(cellType(nc))); % identify target cells
allSlopes = NaN(size(targets));
for nn = 1 : length(targets) 
    prefFreq = Tc.perfFreqPerAmp{targets(nn)}; % preferred freq at each amp
    if sum(prefFreq ~=0) >= 2 % check for at least 2 preferred freq values across amplitude
        trueFreq = NaN(5,1); % yVals (freq)
        for nf = 1 : length(prefFreq)
            if prefFreq(nf) == 0
                continue; %skip this amp
            end
            trueFreq(nf) = freqVal(prefFreq(nf));
        end
        tbl1 = table(ampValAll,trueFreq,'VariableNames',{'x','y'}); % fit mean vals (visualization only)
        f1 = fitlm(tbl1,'y ~ x'); % linear fit
        allSlopes(nn) = f1.Coefficients.Estimate(2); % slope
        % [yPred,yCI] = predict(f1,ampValAll,'Alpha',0.05);
        % plot(ampValAll,yPred,'--','LineWidth',0.5,'Color',colors.(cellType(nc))); % fitted line
    end
end

%histogram(allSlopes,'NumBins',100);
swarmchart(ones(size(allSlopes))+1,allSlopes)
%set(gca,'YScale','log');

%% SELECTIVITY PDF & CDF FROM EFFECT SIZE
% examine cdf of each cell type
figure('theme','light','color',[1 1 1]); hold on;
pThr = 0.05;
for j = 1:length(cellType)
    idx = Tc.identity == cellType(j) & Tc.selectivityPVal < pThr;
    [pdfEffectSize, xValues] = ksdensity(Tc.selectivityEffectSize(idx),'Support',[0 1]);
    plot(xValues, pdfEffectSize, 'DisplayName', char(cellType(j)), 'LineWidth', 2);
end
xlabel('Selectivity Effect Size');
ylabel('Probability Density');
title('PDF of Effect Sizes by Cell Type');
legend('show'); hold off;

% examine cdf of each cell type
figure('theme','light','color',[1 1 1]); hold on;
for j = 1:length(cellType)
    idx = Tc.identity == cellType(j) & Tc.selectivityPVal < pThr;
    [cdfEffectSize, xValuesCDF] = ecdf(Tc.selectivityEffectSize(idx));
    plot(xValuesCDF, cdfEffectSize, 'DisplayName', ['CDF of ' char(cellType(j))], 'LineWidth', 2);
    hold on;
end
xlabel('Selectivity Effect Size');
ylabel('Cumulative probability');
title('CDF of Effect Sizes by Cell Type');
legend('show'); hold off;

%% EVALUATE NEURON TUNINGS (freq & amplitude) WITH PCA (DATASET 2)

tic;
respSpace = zeros(minTrialNum*nFreq*nAmp,size(Tc,1)); % [stim w/repeats (15*5*5) x neurons (11930)]
respSpaceKey = zeros(size(respSpace,1),2); % freq (col.1), amp (col.2)
spx = reshape(1:(minTrialNum*nFreq*nAmp),minTrialNum,nFreq,nAmp);
for nn = 1 : size(Tc,1)
    for nf = 1 : length(freqVal)
        for na = 1 : length(ampValAll)
            idx = find(Tc.frequencies{nn} == freqVal(nf) & Tc.amplitudes{nn} == ampValAll(na));
            locResp = Tc.responses{nn}(idx); % neural responses, this freq & amp
            respSpace(spx(:,na,nf),nn) = locResp(randperm(length(locResp),minTrialNum)); % randomly select 15
            if nn == 1
                respSpaceKey(spx(:,na,nf),1) = repmat(freqVal(nf),minTrialNum,1);
                respSpaceKey(spx(:,na,nf),2) = repmat(ampValAll(na),minTrialNum,1);
            end
        end
    end
end
toc;

% DO MATH
[coeff,score,latent,tsquared] = pca(respSpace(:,is.SOM & is.InterMod));

% PLOT
targetPC = 1;
xValues = reshape(1:(5*9),9,5);
xValues = xValues(1:5,:)';
ampValInv = flip(ampValAll);
figure('theme','light','color',[1 1 1]); hold on;
for n = 1 : 5
    plot(-1,0,'.','Color',colors.Amp(n,:));
end
for na = 1 : length(ampValInv)
    yMeans = zeros(length(freqVal),1);
    for nf = 1 : length(freqVal)
        indx = respSpaceKey(:,1) == freqVal(nf) & respSpaceKey(:,2) == ampValInv(na);
        locResps = score(indx,targetPC);
        plot(repmat(xValues(nf,na),2,1),mkCI(locResps),'-','LineWidth',1.5,'Color',colors.Amp(na,:));
        % % Plot all unique values
        % plot(repmat(xValues(nf,na),minTrialNum,1),locResps,'.','MarkerSize',8,'Color',colors.Amp(na,:));
        yMeans(nf) = mean(locResps);
        plot(xValues(nf,na),yMeans(nf),'.','MarkerSize',12,'Color',colors.Amp(na,:));
    end
    plot(xValues(:,na),yMeans,'-','LineWidth',1,'Color',colors.Amp(na,:))
end
set(gca,'XLim',[0 max(xValues,[],"all")+1],'XTick',xValues(:,3),'XTickLabel',freqVal,...
    'FontName','Arial');
xlabel('Frequency (Hz)','FontName','Arial');
ylabel('\DeltaF/F','FontName','Arial');
title(strcat("Freq/Amp Tuning of PC ",num2str(targetPC)),'FontName','Arial');
legendTxt = {"7.7 \mum","3.5","1.6","0.76","0.33"};
legend(legendTxt);
subtitle("From SOM neurons with Interaction Effect",'FontName','Arial');



%% perform population decoding 

nIter = 200; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
targetAmp = 5; % 0.33,0.76,1.6,3.5,7.7
usePCA = true; % use PCA step prior to model training
kStart = 1; % first PC to include
kComp = 15; % number of PCs to use in model training / testing
useFreqEffNeu = false; % use neurons with significant anova freq effects
useAllNeu = false; % use all available neurons for modeling (no random cell draws)
neuStepSize = 15; % num neurons to add to model each step


if useFreqEffNeu == true
    % PV:446 | SOM:281 | EXC:5596
    numNeuAvail = [sum(is.PV & is.selective); sum(is.SOM & is.selective); sum(is.EXC & is.selective)];
elseif useFreqEffNeu == false
    numNeuAvail = [495; 495; 495];
end

if multiAmp == true
    nCellDraws = 15:5:300;
elseif multiAmp == false
    nCellDraws = 15:neuStepSize:495;
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
        if useFreqEffNeu == true
            if numNeuAvail(ni) < nCellDraws(nd)
                continue
            end
        end

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

nCellDraws = 15:15:495;
figure('Color',[1 1 1],'Theme','Light'); hold on;
% plot fake data for correct legend order
plot([-1 -1],[-1 -1],'k-'); 
plot([-1 -1],[-1 -1],'k--');
fill([-1 -1 -1], [-1 -1 -1],colors.PV,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.SOM,'FaceAlpha',0.5,'EdgeColor','none');
fill([-1 -1 -1], [-1 -1 -1],colors.EXC,'FaceAlpha',0.5,'EdgeColor','none');
for nt = 1 : size(PCTall,4)
    for nc = 1 : length(cellType)
        locMat = squeeze(PCTall(:,nc,:,nt));
        useCols = ~all(locMat == 0,1);
        yMeans = mean(locMat(:,useCols),1).*100;
        % ySems = sem(locMat(:,useCols),1).*100; % Standard error
        yCIs = zeros(sum(useCols),2); % confidence intervals (95%)
        for nd = 1 : sum(useCols)
            yCIs(nd,:) = mkCI(squeeze(PCTall(:,nc,nd,nt))).*100;
        end
        xLooped = [nCellDraws(useCols) flip(nCellDraws(useCols))]; 
        % yLooped = [yMeans+ySems flip(yMeans-ySems)]; % when using SEM
        yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
        fill(xLooped,yLooped,colors.(cellType(nc)),'FaceAlpha',0.5,'EdgeColor','none');
        switch nt 
            case 1
                plot(nCellDraws(useCols),yMeans,'-','LineWidth',1,'Color',colors.(cellType(nc)));
            case 2
                plot(nCellDraws(useCols),yMeans,'--','LineWidth',1,'Color',colors.(cellType(nc)));
        end
    end
end
% plot([15 500],[20 20],'r--'); % chance performance level 
% set(gca,'XLim',[15 500],'XTick',[15:50:500 500],'YLim',[30 100],'YTick',30:10:100);
set(gca,'XLim',[15 500],'XTick',[15:50:500 500],'YLim',[25 70],'YTick',20:10:100);
ylabel('model accuracy (%)','FontName','Arial');
xlabel('Neurons used','FontName','Arial');
title('Frequency decoding performance by cell type','FontName','Arial');
subtitle('From neurons with frequency selectivity','FontName','Arial');
legend({"Freq. Mod.","All neurons","PV","SOM","EXC"},'Location','NorthWest');


%% PLOT ACCURACY TOGETHER (old way?)
% figure('Color',[1 1 1],'Theme','Light');
% v = violinplot(PCTall); hold on;
% v(1).FaceColor = colors.PV; v(2).FaceColor = colors.SOM; v(3).FaceColor = colors.EXC;
% v(4).FaceColor = colors.PV; v(5).FaceColor = colors.SOM; v(6).FaceColor = colors.EXC;
% v(7).FaceColor = colors.PV; v(8).FaceColor = colors.SOM; v(9).FaceColor = colors.EXC;
% v(10).FaceColor = colors.PV; v(11).FaceColor = colors.SOM; v(12).FaceColor = colors.EXC;
% if multiAmp == true
%     set(gca,'XTick',["2","5","8"],'XTickLabels',string(nCellDraws),'FontName','Arial');
% elseif multiAmp == false
%     set(gca,'XTick',["2","5","8","11"],'XTickLabels',string(nCellDraws),'FontName','Arial');
% end
% if usePCA == true
%     title("Frequency decoding through LDA",'FontName','Arial');
% elseif usePCA == false
%     title("LDA classification without PCA",'FontName','Arial');
% end
% 
% if useFreqEffNeu == true
%     if multiAmp == true
%         subtitle({strcat("Model iterations: ",num2str(nIter),", training ratio: ",num2str(trainPortion)),...
%             strcat("Amplitude: ",num2str(ampValAll(targetAmp)),"um, of Freq Effect Neu.")},'FontName','Arial');
%     elseif multiAmp == false
%         subtitle(strcat("Model iterations: ",num2str(nIter),", training ratio: ",num2str(trainPortion)),'FontName','Arial')
%     end
% elseif useFreqEffNeu == false
%     if multiAmp == true
%         subtitle({strcat("Model iterations: ",num2str(nIter),", training ratio: ",num2str(trainPortion)),...
%             strcat("Amplitude: ",num2str(ampValAll(targetAmp)),"um")},'FontName','Arial');
%     elseif multiAmp == false
%         % no subtitle
%     end
% end
% xlabel('neurons per model','FontName','Arial')
% ylabel('decoding accuracy','FontName','Arial');
% set(gca,'YLim',[-0 1.1]);
% plot([3.5 3.5],[-0.1 1.1],'k--');
% plot([6.5 6.5],[-0.1 1.1],'k--');
% plot([9.5 9.5],[-0.1 1.1],'k--');
% 
% if includeErrBar == true
%     semsPV = [sem(PCTall(:,1)), sem(PCTall(:,4)), sem(PCTall(:,7)), sem(PCTall(:,10))];
%     errorbar([1 4 7 10],mean(PCTall(:,[1 4 7 10]),1),semsPV,"LineStyle","none",'Color',colors.PV,'CapSize',10);
% 
%     semsSOM = [sem(PCTall(:,2)), sem(PCTall(:,5)), sem(PCTall(:,8)), sem(PCTall(:,11))];
%     errorbar([2 5 8 11],mean(PCTall(:,[2 5 8 11]),1),semsSOM,"LineStyle","none",'Color',colors.SOM,'CapSize',10);
% 
%     semsEXC = [sem(PCTall(:,3)), sem(PCTall(:,6)), sem(PCTall(:,9)), sem(PCTall(:,12))];
%     errorbar([3 6 9 12],mean(PCTall(:,[3 6 9 12]),1),semsEXC,"LineStyle","none",'Color',colors.EXC,'CapSize',10);
% end
% 
% legend({"PV","SOM","EXC"},'Location','SouthEast');
% 
% % % SAVE
% % saveas(gcf,'dataset1_LDA_performance','svg')
% 
% % PLOT CONFUSION MATRICIES
% for nf = 1 : length(cellType)
%     locPreds = squeeze(predFreqAll(:,nf,:));
%     locTrues = squeeze(testFreqAll(:,nf,:));
%     figure('Theme','light','Position',[100 100 400 400]);
%     confusionchart(locPreds(:),locTrues(:),'DiagonalColor',colors.(cellType(nf)),'OffDiagonalColor',[0.15 0.15 0.15]);
%     title(strcat("Confusion matrix: ",cellType(nf)));
%     set(gca,'FontName','Arial');
%     saveas(gcf,strcat("dataset1_LDA_confuMat_",cellType(nf)),'svg')
% end


%% PLOT PC TUNING CURVES FROM MODELS
targetCellType = 1; %(PV, SOM, EXC)
targetDrawBin = 12;
RGB = mkcolors([1 0 0],[0 0 1],10);
figure('Theme','light','color',[1 1 1]);
tiledlayout(1,5);
for nr = 1 : 5 %num reps
mdl = MDLall{nr,targetCellType,targetDrawBin}; % nIter (200), nCellTypes (3), nCellDraws (4)
PCresps = zeros(20,5,15);
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
if nr == 1
legend({'PC1','','','','','','PC2','','','','','','PC3','','','','','','PC4',...
    '','','','','','PC5','','','','','','PC6','','','','','','PC7'},'Location','southeast');
end
hold off;
end


%% PLOT PC TUNING CURVES AS HEAT MAPS
targetDrawBin = 18;
targetRep = 1;
numPCs = 15;
RGB = mkcolors([1 0 0],[0 0 1],numPCs);
figure('Theme','light','color',[1 1 1]);
tiledlayout(1,length(cellType));
for nc = 1 : length(cellType)
    mdl = MDLall{targetRep,nc,targetDrawBin,1}; % nIter (200), nCellTypes (3), nCellDraws (4)
    PCresps = zeros(20,5,numPCs);
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


%% Evaluate MDL PC selectivity (new way)
% top 5 PCs, each cell type compared
figure('theme','light','color',[1 1 1]); hold on;
numPCEval = 6; %num PCs to evaluate
targetDrawBin = 1; % TOP ONLY
PCresps = zeros(minTrialNum*trainPortion,nFreq,numPCEval); % 20x5x5
locSel = zeros(numReps,1);
xp = 0.1;
for np = 1 : numPCEval
    for nc = 1 : length(cellType)
        for nr = 1 : numReps
            mdl = MDLall{nr,nc,targetDrawBin};
            for nf = 1 : nFreq
                PCresps(:,nf,np) = mdl.X(mdl.Y == freqVal(nf),np);
            end
            locSel(nr) = quickSelectivity(mean(squeeze(PCresps(:,:,np)),1));
        end
        switch nc
            case 1 % PV
                plot(np-xp,mean(locSel),'.','MarkerSize',8,'Color',colors.PV);
                plot([np np]-xp,[mean(locSel)-sem(locSel) mean(locSel)+sem(locSel)],'-','Color',colors.PV);
            case 2 % SOM
                plot(np,mean(locSel),'.','MarkerSize',8,'Color',colors.SOM);
                plot([np np],[mean(locSel)-sem(locSel) mean(locSel)+sem(locSel)],'-','Color',colors.SOM);
            case 3 % EXC
                plot(np+xp,mean(locSel),'.','MarkerSize',8,'Color',colors.EXC);
                plot([np np]+xp,[mean(locSel)-sem(locSel) mean(locSel)+sem(locSel)],'-','Color',colors.EXC);
        end
    end
end
set(gca,'XTick',1:numPCEval);
xlabel("PC number",'FontName','Arial');
ylabel("Selectivity Strength",'FontName','Arial');
title("Frequency selectivity strength of top PCs across models",'FontName','Arial');
legend({"PV","","SOM","","EXC",""},'FontName','Arial');




%% Evaluate MDL PC selectivity, one subplot per cell type (old way)

for nc = 1 : 3 % LOOP: CELL TYPE (3)
    nexttile(nc); hold on;
    targetCellType = nc;
    pcEvalN = 5; %number of PCs for evaluation
    numReps = 200;
    numDrawBins = 1;
    % numTotalBins = 1;
    % multiDrawBins = round(linspace(1,numTotalBins,numDrawBins));
    % RGB = mkcolors([1 0 0],[0 0 1],numDrawBins);
    xp = 0.1; % x-padding
    for nb = 1 : numDrawBins %LOOP: DRAW BINS (1)
        locSel = zeros(numReps,pcEvalN);
        targetDrawBin = nb;
        xpl = xp+(nb-(round(numDrawBins/2)+1))*xp;
        for nr = 1 : numReps % LOOP: NUM REPS (200)
            mdl = MDLall{nr,targetCellType,targetDrawBin};
            PCresps = zeros(minTrialNum*trainPortion,nFreq,pcEvalN); % 20x5x15
            for np = 1 : pcEvalN 
                for nf = 1 : nFreq
                    PCresps(:,nf,np) = mdl.X(mdl.Y == freqVal(nf),np);
                end
                locSel(nr,np) = quickSelectivity(mean(squeeze(PCresps(:,:,np)),1));
            end
        end
        locSems = sem(locSel);
        locMeans = mean(locSel);
        for np = 1 : pcEvalN
            yVals = [locMeans(np)+locSems(np) locMeans(np)-locSems(np)];
            plot([np np]+xpl,yVals,'-','Color',RGB(nb,:));
        end
        plot((1:pcEvalN)+xpl,locMeans,'.','Color',RGB(nb,:),'MarkerSize',10);
    end
    hold off;
    if nc == 3
        xlabel("PC number",'FontName','Arial');
        set(gca,'XTick',1:15);
    else
        set(gca,'XTick',1:15,'XTickLabel',[]);
    end
    if nc == 1 
        title("Selectivity strength of PCs as models increase",'FontName','Arial');
    end
    ylabel("selectivity strength",'FontName','Arial');
    set(gca,'YLim',[0.25 0.4]);
    tTxt = strcat(cellType(targetCellType)," neuron models");
    subtitle(tTxt,'FontName','Arial');
    hold off;
end


%% Evaluate MDL PC response magnitude
figure('theme','light','Color',[1 1 1]); 
tiledlayout(3,1);
for nc = 1 : 3
nexttile(nc); hold on;
targetCellType = nc;
pcEvalN = 15; %number of PCs for evaluation
numReps = 200;
numDrawBins = 5;
numTotalBins = 32;
multiDrawBins = round(linspace(1,numTotalBins,numDrawBins));
RGB = mkcolors([1 0 0],[0 0 1],numDrawBins);
xp = 0.1; % x-padding
for nb = 1 : numDrawBins
    % locSel = zeros(numReps,pcEvalN);
    locMag = zeros(numReps,pcEvalN);
    targetDrawBin = nb;
    xpl = xp+(nb-(round(numDrawBins/2)+1))*xp;
    for nr = 1 : numReps
        mdl = MDLall{nr,targetCellType,targetDrawBin};
        PCresps = zeros(minTrialNum*trainPortion,nFreq,pcEvalN); % 20x5x15
        for np = 1 : pcEvalN
            for nf = 1 : nFreq
                PCresps(:,nf,np) = mdl.X(mdl.Y == freqVal(nf),np);
            end
            % locSel(nr,np) = quickSelectivity(mean(squeeze(PCresps(:,:,np)),1));
            avgResps = mean(squeeze(PCresps(:,:,np)),1);
            locMag(nr,np) = abs(diff([max(avgResps) min(avgResps)]));
        end
    end
    locSems = sem(locMag);
    locMeans = mean(locMag);
    for np = 1 : pcEvalN
        yVals = [locMeans(np)+locSems(np) locMeans(np)-locSems(np)];
        plot([np np]+xpl,yVals,'-','Color',RGB(nb,:));
    end
    plot((1:pcEvalN)+xpl,locMeans,'.','Color',RGB(nb,:),'MarkerSize',10);
end
hold off;
if nc == 3
    xlabel("PC number",'FontName','Arial');
    set(gca,'XTick',1:15);
else
    set(gca,'XTick',1:15,'XTickLabel',[]);
end
if nc == 1 
    title("Frequency response difference of PCs as models increase",'FontName','Arial');
end
ylabel("max response difference",'FontName','Arial');
%set(gca,'YLim',[0.25 0.4]);
tTxt = strcat(cellType(targetCellType)," neuron models");
subtitle(tTxt,'FontName','Arial');
hold off;
end



%% Evaluate MDL PC frequency distribution
figure('theme','light','Color',[1 1 1]); 
tiledlayout(3,15);
xVals = reshape(1:(7*15),[7 15])';
xVals()
for nc = 1 : 3
targetCellType = nc;
pcEvalN = 15; %number of PCs for evaluation
numReps = 200;
numDrawBins = 5;
numTotalBins = 32;
multiDrawBins = round(linspace(1,numTotalBins,numDrawBins));
RGB = mkcolors([1 0 0],[0 0 1],numDrawBins);
xp = 0.1; % x-padding
allPCT = zeros(15,5,5); % PC num, freq, numDrawBins
for nb = 1 : numDrawBins
    % locSel = zeros(numReps,pcEvalN);
    locFreq = zeros(numReps,pcEvalN);
    targetDrawBin = nb;
    xpl = xp+(nb-(round(numDrawBins/2)+1))*xp;
    for nr = 1 : numReps
        mdl = MDLall{nr,targetCellType,targetDrawBin};
        PCresps = zeros(minTrialNum*trainPortion,nFreq,pcEvalN); % 20x5x15
        for np = 1 : pcEvalN
            for nf = 1 : nFreq
                PCresps(:,nf,np) = mdl.X(mdl.Y == freqVal(nf),np);
            end
            % locSel(nr,np) = quickSelectivity(mean(squeeze(PCresps(:,:,np)),1));
            avgResps = mean(squeeze(PCresps(:,:,np)),1);
            [~,locFreq(nr,np)] = max(abs(avgResps));
        end
    end
    % locSems = sem(locMag);
    % locMeans = mean(locMag);
    % calculate percentage of 1-5 in each column of locFreq
    allPCT(:,:,nb) = 100 * squeeze(sum(locFreq == reshape(1:5,1,1,[]), 1)) / size(locFreq,1);
    % for np = 1 : pcEvalN
    %     yVals = [locMeans(np)+locSems(np) locMeans(np)-locSems(np)];
    %     plot([np np]+xpl,yVals,'-','Color',RGB(nb,:));
    % end
    % plot((1:pcEvalN)+xpl,locMeans,'.','Color',RGB(nb,:),'MarkerSize',10);
end

for np = 1 : pcEvalN
    nexttile; 
    bar(squeeze(allPCT(np,:,:))','stacked');
    if nc == 1 
        subtitle(strcat("PC ",num2str(np)),'FontName','Arial');
        set(gca,'XTickLabel',[]);
        if np == 8
            title("Frequency Preference of PCs across all reps",'FontName','Arial');
        end
        if np == 15
            legend({"100Hz","300","500","700","900","1100"});
        end
    end
    if nc == 2 & np == 1
        ylabel("percentage (%)",'FontName','Arial');
    end
    if nc == 3 & np == 8
        xlabel("Neuron bin number",'FontName','Arial');
    end
    hold off
end



end




%%  ONE TIME
% perform population decoding


nIter = 100; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
targetAmp = 5; % 0.33,0.76,1.6,3.5,7.7
usePCA = true; % use PCA step prior to model training
kComp = 15; %number of PCs to use in model training / testing
useFreqEffNeu = true; % use neurons with significant anova freq effects
useAllNeu = true; % use all available neurons for modeling
includeErrBar = true; %include error bars in plots
cellType = ["PV"; "SOM"; "EXC"]; % order for plotting




tic;
for ni = 1:length(cellType)

    disp(cellType(ni))

    % populate arrays
    nCellDraw = nCellDraws(nd);
    trainTrial = round(minTrialNum*trainPortion);
    testTrial = round(minTrialNum*(1-trainPortion));
    trainFreq = repmat(freqVal,trainTrial,1);
    testFreq = repmat(freqVal,testTrial,1);
    PerfCellType = zeros(nIter,nCellType);
    predFreqAll = zeros(nIter,nCellType,nFreq*testTrial);
    testFreqAll = zeros(nIter,nCellType,nFreq*testTrial);
    LR_rsqs = zeros(nIter,nCellType);

    % select cell type
    if useFreqEffNeu == true
        if multiAmp == true
            idx = find(Tc.identity == cellType(ni) & is.InterMod);
        elseif multiAmp == false
            idx = find(Tc.identity == cellType(ni) & Tc.selectivityPVal < 0.05);
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
            % kComp = 20;
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
        PerfCellType(nj,ni) = acc;

        figure;
        confusionchart(testFreqAll(:),predFreqAll(:))

    end
end

