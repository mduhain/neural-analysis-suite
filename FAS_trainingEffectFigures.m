% FAS_trainingEffectFigures.m
%
% From focusedAnalysisScript loading "dataset2_updated_20260202.mat'
% 
% 2026-02-02 mduhain: Initial practice run with traind AmpFreqMod data
%
%
%

highFreqMice = [3 4 8 12]; % from mouseNums
lowFreqMice = [1 7 2 11];



%% calculate preferred frequency per amplitude, selectivity & sel. metric (strength)

pValThresh = 0.05;
for nn = 1 : size(Tc,1) % Loop neurons
    % preallocate variables for each neuron
    locResps = Tc.responses{nn}; % response values per trial
    locFreqs = Tc.frequencies{nn}; % frequency values per trial
    locAmps = Tc.amplitudes{nn}; % amplitude values per trial
    Tc.meanRespAF{nn} = zeros(length(ampValAll),length(freqVal)); % amps (5) x freqs (5)
    Tc.medianRespAF{nn} = zeros(length(ampValAll),length(freqVal)); % amps (5) x freqs (5)
    Tc.perfFreqPerAmp{nn} = zeros(length(ampValAll),1); % non-selective (0) or freq. sel (1-5)
    Tc.selMetricPerAmp{nn} = zeros(length(ampValAll),1); % selectivity strength per amp (5x1)
    for na = 1 : length(ampValAll) % loop amplitudes
        allResps = cell(length(freqVal),1); % responses per freq (5x1 cell)
        allFreqs = cell(length(freqVal),1); % freq labels (5x1 cell)
        for nf = 1 : length(freqVal) % loop frequencies
             allResps{nf} = locResps(locFreqs == freqVal(nf) & locAmps == ampValAll(na));
             allFreqs{nf} = repmat(nf,length(allResps{nf}),1);
        end
        Tc.meanRespAF{nn}(na,:) = cellfun(@mean,allResps)'; % extract mean response per freq
        Tc.medianRespAF{nn}(na,:) = cellfun(@median,allResps)'; % median response per freq
        [p,~] = anova1(cell2mat(allResps),cell2mat(allFreqs),'off'); % one-way anova for selectivity
        if p < pValThresh % Neuron has freq. selectivity
            [~,Tc.perfFreqPerAmp{nn}(na)] = max(Tc.meanRespAF{nn}(na,:)); % save pref. freq.
            Tc.selMetricPerAmp{nn}(na) = quickSelectivity(Tc.meanRespAF{nn}(na,:)); % selectivity strength
        end
    end
end

clear locResps locFreqs locAmps allResps allFreqs p tbl nn na nf
% save('dataset2_forFigures.mat','-v7.3');

%% NEURON COUNTS, BY CELL TYPE, AND BY SELECTIVITY TYPE

% SAVED IN ABOVE (2025-11-17)^^
num = struct(); % structure to hold neuron counts from all type combinations

% PV NEURONS COUNTS
num.PV = sum(is.PV);
num.unmodPV = sum(is.PV & is.nonMod); % unmodulated
num.freqPV = sum(is.PV & is.FreqMod & ~is.AmpMod & ~is.InterMod); % freq effect alone
num.ampPV = sum(is.PV & is.AmpMod & ~is.FreqMod & ~is.InterMod); % amp effect alone
num.dualPV = sum(is.PV & is.DualMod & ~is.InterMod); % dual modulated (2 main effects)
num.interPV = sum(is.PV & is.InterMod); % interaction effect

% SOM NEURON COUNTS
num.SOM = sum(is.SOM);
num.unmodSOM = sum(is.SOM & is.nonMod); % unmodulated
num.freqSOM = sum(is.SOM & is.FreqMod & ~is.AmpMod & ~is.InterMod); % freq effect alone
num.ampSOM = sum(is.SOM & is.AmpMod & ~is.FreqMod & ~is.InterMod); % amp effect alone
num.dualSOM = sum(is.SOM & is.DualMod & ~is.InterMod); % dual modulated (2 main effects)
num.interSOM = sum(is.SOM & is.InterMod); % interaction effect

% EXC NEURON COUNTS
num.EXC = sum(is.EXC);
num.unmodEXC = sum(is.EXC & is.nonMod); % unmodulated
num.freqEXC = sum(is.EXC & is.FreqMod & ~is.AmpMod & ~is.InterMod); % freq effect alone
num.ampEXC = sum(is.EXC & is.AmpMod & ~is.FreqMod & ~is.InterMod); % amp effect alone
num.dualEXC = sum(is.EXC & is.DualMod & ~is.InterMod); % dual modulated (2 main effects)
num.interEXC = sum(is.EXC& is.InterMod); % interaction effect

%% PLOT % RESPONSIVITY PRE VS POST
figure('theme','light','color',[1 1 1]); hold on;
respPerc = zeros(length(cellType),2);
xVals = [1 2 4 5 7 8];
for i = 1:length(cellType)
    % RESPONSIVE UNITS CALC
    indxPre = Tc.identity == cellType(i) & Tc.isTrained == false; % target neurons (pre-training)
    indxPost = Tc.identity == cellType(i) & Tc.isTrained == true; % target neurons (post-training)
    numRespsPre = sum(Tc.responsivityPVal(indxPre)<0.05); % number of responsive neurons (pre-training)
    numRespsPost = sum(Tc.responsivityPVal(indxPost)<0.05); % number of responsive neurons (post-training)
    respPerc(i,1) = numRespsPre/sum(indxPre); % responsive neuron percentage (pre-training)
    respPerc(i,2) = numRespsPost/sum(indxPost); % responsive neuron percentage (post-training)
    respSEPre = mannyPropSTD(numRespsPre,sum(indxPre)); % responsive neuron estimated standard error
    respSEPost = mannyPropSTD(numRespsPost,sum(indxPost)); % responsive neuron estimated standard error
    respCIPre = respSEPre * 1.96; %SE to 95% CI
    respCIPost = respSEPost * 1.96; %SE to 95% CI
    
    % RESPOINSIVE UNITS PLOT (pre-training)
    bar(xVals(i*2-1),respPerc(i,1)*100,'FaceColor',colors.(cellType(i))); 
    text(xVals(i*2-1),((respPerc(i,1)*100)+(respCIPre*100)+0.02),strcat("n=",num2str(numRespsPre)),...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Arial','FontSize',9);
    errorbar(xVals(i*2-1),respPerc(i,1)*100,respCIPre*100,Color=[0,0,0],LineWidth=1.4);
    text(xVals(i*2-1),71,'Before','HorizontalAlignment','left','VerticalAlignment','middle',...
        'Rotation',90,'FontName','Arial','FontSize',10,'Color',[1 1 1]);

    % RESPOINSIVE UNITS PLOT (post-training)
    bar(xVals(i*2),respPerc(i,2)*100,'FaceColor',colors.(cellType(i))); 
    text(xVals(i*2),((respPerc(i,2)*100)+(respCIPost*100)+0.02),strcat("n=",num2str(numRespsPost)),...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Arial','FontSize',9)
    errorbar(xVals(i*2),respPerc(i,2)*100,respCIPost*100,Color=[0,0,0],LineWidth=1.4);
    text(xVals(i*2),71,'After','HorizontalAlignment','left','VerticalAlignment','middle',...
        'Rotation',90,'FontName','Arial','FontSize',10,'Color',[1 1 1]);

end
set(gca,'XTick',[1.5 4.5 7.5],'XTickLabels',cellType,'ylim',[70 102],'FontName','Arial','Box','off');
title('Training effect on stimulus responsivity','FontName','Arial');
ylabel('Percentage of neurons per group','FontName','Arial');


%% PLOT % FREQUENCY SELECTIVITY (F)
figure('theme','light','color',[1 1 1]); hold on;
selPerc = zeros(length(cellType),2);
xVals = [1 2 4 5 7 8];
for i = 1:length(cellType)
    % SELECTIVE UNITS CALC
    indxPre = Tc.identity == cellType(i) & Tc.isTrained == false; % target neurons
    indxPost = Tc.identity == cellType(i) & Tc.isTrained == true; % target neurons
    numRespsPre = sum(Tc.responsivityPVal(indxPre)<0.05); % number of responsive neurons
    numRespsPost = sum(Tc.responsivityPVal(indxPost)<0.05); % number of responsive neurons
    numSelsPre = sum(Tc.anovanFreqP(indxPre)<0.05); % number of frequency-effect (F) neurons
    numSelsPost = sum(Tc.anovanFreqP(indxPost)<0.05); % number of frequency-effect (F) neurons
    selPerc(i,1) = numSelsPre/numRespsPre;
    selPerc(i,2) = numSelsPost/numRespsPost;
    selSEPre = mannyPropSTD(numSelsPre,numRespsPre); % selective neuron estimated standard error
    selSEPost = mannyPropSTD(numSelsPost,numRespsPost); % selective neuron estimated standard error
    selCIPre = selSEPre * 1.96;%SE to 95% CI
    selCIPost = selSEPost * 1.96;%SE to 95% CI

    % SELECTIVE UNITS PLOT (pre-training)
    bar(xVals(i*2-1),selPerc(i,1)*100,'FaceColor',colors.(cellType(i))); 
    text(xVals(i*2-1),((selPerc(i,1)*100)+(selCIPre*100)+0.02),strcat("n=",num2str(numSelsPre)),'HorizontalAlignment',...
        'center','VerticalAlignment','bottom','FontName','Arial','FontSize',9)
    errorbar(xVals(i*2-1),selPerc(i,1)*100,selCIPre*100,Color=[0,0,0],LineWidth=1.4);
    text(xVals(i*2-1),15,'Before','HorizontalAlignment','left','VerticalAlignment','middle',...
        'Rotation',90,'FontName','Arial','FontSize',10,'Color',[1 1 1]);

    % SELECTIVE UNITS PLOT (post-training)
    bar(xVals(i*2),selPerc(i,2)*100,'FaceColor',colors.(cellType(i))); 
    text(xVals(i*2),((selPerc(i,2)*100)+(selCIPost*100)+0.02),strcat("n=",num2str(numSelsPost)),'HorizontalAlignment',...
        'center','VerticalAlignment','bottom','FontName','Arial','FontSize',9)
    errorbar(xVals(i*2),selPerc(i,2)*100,selCIPost*100,Color=[0,0,0],LineWidth=1.4);
    text(xVals(i*2),15,'After','HorizontalAlignment','left','VerticalAlignment','middle',...
        'Rotation',90,'FontName','Arial','FontSize',10,'Color',[1 1 1]);

end
set(gca,'Xtick',[1.5 4.5 7.5],'XTickLabels',cellType,'ylim',[14 83],'FontName','Arial','Box','off');
title('Training effect on feature interaction (F)','FontName','Arial');
ylabel('Percentage of responsive neurons','FontName','Arial');

%% PLOT % FEATURE INTERACTION (F*A)
figure('theme','light','color',[1 1 1]); hold on;
selPerc = zeros(length(cellType),2);
xVals = [1 2 4 5 7 8];
for i = 1:length(cellType)
    % SELECTIVE UNITS CALC
    indxPre = Tc.identity == cellType(i) & Tc.isTrained == false; % target neurons
    indxPost = Tc.identity == cellType(i) & Tc.isTrained == true; % target neurons
    numRespsPre = sum(Tc.responsivityPVal(indxPre)<0.05); % number of responsive neurons
    numRespsPost = sum(Tc.responsivityPVal(indxPost)<0.05); % number of responsive neurons
    numSelsPre = sum(Tc.anovanMixP(indxPre)<0.05); % number of interaction-effect (F*A) neurons
    numSelsPost = sum(Tc.anovanMixP(indxPost)<0.05); % number of interaction-effect (F*A) neurons
    selPerc(i,1) = numSelsPre/numRespsPre;
    selPerc(i,2) = numSelsPost/numRespsPost;
    selSEPre = mannyPropSTD(numSelsPre,numRespsPre); % selective neuron estimated standard error
    selSEPost = mannyPropSTD(numSelsPost,numRespsPost); % selective neuron estimated standard error
    selCIPre = selSEPre * 1.96;%SE to 95% CI
    selCIPost = selSEPost * 1.96;%SE to 95% CI

    % SELECTIVE UNITS PLOT (pre-training)
    bar(xVals(i*2-1),selPerc(i,1)*100,'FaceColor',colors.(cellType(i))); 
    text(xVals(i*2-1),((selPerc(i,1)*100)+(selCIPre*100)+0.02),strcat("n=",num2str(numSelsPre)),'HorizontalAlignment',...
        'center','VerticalAlignment','bottom','FontName','Arial','FontSize',9)
    errorbar(xVals(i*2-1),selPerc(i,1)*100,selCIPre*100,Color=[0,0,0],LineWidth=1.4);
    text(xVals(i*2-1),1,'Before','HorizontalAlignment','left','VerticalAlignment','middle',...
        'Rotation',90,'FontName','Arial','FontSize',10,'Color',[1 1 1]);

    % SELECTIVE UNITS PLOT (post-training)
    bar(xVals(i*2),selPerc(i,2)*100,'FaceColor',colors.(cellType(i))); 
    text(xVals(i*2),((selPerc(i,2)*100)+(selCIPost*100)+0.02),strcat("n=",num2str(numSelsPost)),'HorizontalAlignment',...
        'center','VerticalAlignment','bottom','FontName','Arial','FontSize',9)
    errorbar(xVals(i*2),selPerc(i,2)*100,selCIPost*100,Color=[0,0,0],LineWidth=1.4);
    text(xVals(i*2),1,'After','HorizontalAlignment','left','VerticalAlignment','middle',...
        'Rotation',90,'FontName','Arial','FontSize',10,'Color',[1 1 1]);

end
set(gca,'Xtick',[1.5 4.5 7.5],'XTickLabels',cellType,'ylim',[0 40],'FontName','Arial','Box','off');
title('Training effect on feature interaction (F*A)','FontName','Arial');
ylabel('Percentage of responsive neurons','FontName','Arial');


%% Percentage of freq preferring neurons per amp
% as percentage of all neurons with selectivity

selType = ["FreqMod","DualMod","InterMod"];
% selType = ["FreqMod"];
% ns = 1;
selTypeName = ["F","F+A","F*A"];
% selTypeName = "F";
perfFreqGrid = cat(2,Tc.perfFreqPerAmp{:}); % 0 == NonSelective

nc = 1; %target cell type (PV, SOM, EXC)

xSp = 2; % X spacing for plotting
nBins = 2; %numel(cellType); 
nTypes = length(ampValAll); 
d1 = nTypes + xSp; 
xMat = (1:d1*nBins);
xMat = reshape(xMat, d1, nBins).'; 
xMat = xMat(:, 1:nTypes); 
xMat = xMat';

for ns = 1 : length(selType)

    figure('theme','light','color',[1 1 1]); hold on;
    allYs = [];
    
    % pre training
    targets = is.(selType(ns)) & is.(cellType(nc)) & ~is.Trained & any(is.mouseNum(:,highFreqMice),2);
    for na = 1 : length(ampValAll)
        locVal = zeros(length(freqVal),1);
        for nf = 1 : length(freqVal)
            locVal(nf) = sum(perfFreqGrid(na,targets) == nf);
            % locVal(nf) = sum(perfFreqGrid(na,is.(selType(ns))) == nf);
        end
        allYs = cat(1,allYs,(locVal./sum(locVal))');
    end
    text(xMat(3,1),102,'Before','FontName','Arial','FontSize',11,'HorizontalAlignment','center');
    
    % post training
    targets = is.(selType(ns)) & is.(cellType(nc)) & is.Trained & any(is.mouseNum(:,highFreqMice),2);
    for na = 1 : length(ampValAll)
        locVal = zeros(length(freqVal),1);
        for nf = 1 : length(freqVal)
            locVal(nf) = sum( perfFreqGrid(na,targets) == nf );
            % locVal(nf) = sum(perfFreqGrid(na,is.(selType(ns))) == nf);
        end
        allYs = cat(1,allYs,(locVal./sum(locVal))');
    end
    
    text(xMat(3,2),102,'After','FontName','Arial','FontSize',11,'HorizontalAlignment','center');
    
    b = bar(xMat(:),allYs.*100,'stacked');
    cmap = colormap("parula");
    faceColors = cmap(round(linspace(1,256,5)),:); % matlab default
    % Set the face colors for each bar
    for k = 1:length(b)
        b(k).FaceColor = 'flat'; % Use flat color for individual bars
        b(k).CData = faceColors(k, :); % Assign the specified color
    end
    
    set(gca,'XTick',xMat(:),'XTickLabels',repmat(ampValAll,3,1),'XTickLabelRotation',45,...
        'YLim',[0 105],'YTickLabels',0:10:100,'XLim',[0 max(xMat(:))+xSp*2],'FontName','Arial');
    legendTxt = string(freqVal);
    legendTxt(1) = strcat(legendTxt(1)," Hz");
    legend(legendTxt,'box','off','Location','east');
    ylabel('Percentage of selective neurons','FontName','Arial');
    xlabel("Amplitude (\mum)",'FontName','Arial');
    % title(strcat(cellType(nc)," (",selTypeName(ns),") preferred frequency by amplitude"),'FontName','Arial');
    title("Effect of high-frequency training",'FontName','Arial');
    subtitle(strcat("Preferred frequency of (",selTypeName(ns),") ",cellType(nc)," neurons"),'FontName','Arial');
end





%% Number of selective neurons by amplitude (MEAN POINTS + SE)
% as percentage of all neurons with selectivity

num.fePV = sum(is.PV & (is.FreqMod | is.DualMod | is.InterMod));
num.feSOM = sum(is.SOM & (is.FreqMod | is.DualMod | is.InterMod));
num.feEXC = sum(is.EXC & (is.FreqMod | is.DualMod | is.InterMod));
selMetricGrid = cat(2,Tc.selMetricPerAmp{:}); % 0 == NonSelective

xSp = 0.1; % X spacing for plotting
xMat = repmat((1:5)',1,length(cellType))+repmat([xSp*-1 0 xSp],length(ampValAll),1);
yMat = zeros(size(xMat));

figure('theme','light','color',[1 1 1]); hold on;

% blank data for correct legend order
for nc = 1 : length(cellType)
    plot(-1,-1,'.','MarkerSize',10,'Color',colors.(cellType(nc)));
end
plot(-1,-1,'--','LineWidth',0.5,'Color',[0 0 0]);

for nc = 1 : length(cellType) 
    for na = 1 : length(ampValAll)
        locVal = sum(selMetricGrid(na,is.(cellType(nc))) ~= 0);
        yVal = locVal/num.('fe'+cellType(nc)) * 100; % percentage
        yMat(na,nc) = yVal;
        plot(xMat(na,nc),yVal,'.','MarkerSize',10,'Color',colors.(cellType(nc)));
        locSE = mannyPropSTD(locVal,num.('fe'+cellType(nc))) * 100;
        locCI = locSE*1.96;
        errorbar(xMat(na,nc),yVal,locCI,'vertical','Color',colors.(cellType(nc)),'LineWidth',1);
        % text(xMat(na,nc),20,strcat(" n = ",num2str(locVal)),'HorizontalAlignment','left',...
        %     'VerticalAlignment','middle','FontName','Arial','FontSize',9,'Rotation',90,'Color',[0 0 0]);
    end
    tbl = table(xMat(:,nc),yMat(:,nc),'VariableNames',{'x','y'});
    f1 = fitlm(tbl,'y ~ x'); %linear fit to plotted values (visualiation only)
    [yPred,yCI] = predict(f1,xMat(:,nc),'Alpha',0.05);
    plot(xMat(:,nc),yPred,'--','LineWidth',0.5,'Color',colors.(cellType(nc))); % fitted line
end


set(gca,'XTick',xMat(:,2),'XLim',[xMat(1)-(xSp*2) xMat(end)+(xSp*2)],'XTickLabels',ampValAll,...
    'YLim',[20 70],'FontName','Arial');
ylabel({'Percentage of frequency selective neurons','(F, F+A, or F*A)'},'FontName','Arial');
xlabel('Amplitude (\mum)','FontName','Arial');
title({'High amplitudes increase frequency','selectivity abundance'},'FontName','Arial');
% subtitle('From neurons any frequency modulation within each cell type','FontName','Arial');
legend({'PV','SOM','EXC'},'Location','northwest','Box','off');


p_values = zeros(length(cellType),1);
% PERMUTATION TESTS
for nc = 1 : length(cellType)
    b_obs = polyfit(ampValAll, yMat(:,nc), 1); % fit amplitude x percentage slope
    b_obs = b_obs(1); % slope
    P = perms(yMat(:,nc)); % permutations (shuffle 5 groups order, 120 permutations)
    nPerm = size(P,1); % 120
    b_perm = zeros(nPerm,1); %
    for np = 1:nPerm
        b = polyfit(ampValAll,P(np,:)',1);
        b_perm(np) = b(1);
    end
    p_values(nc) = mean(abs(b_perm) >= abs(b_obs)); % p-value (two sided test)
end




































