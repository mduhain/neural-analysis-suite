% FAS_selTypePerf_allCells_ampAndFreq.m
%
% FOr figure 7 of manuscript:
% Decode freq then amp sequentially from the same set of neurons
% Make comparisons between how each neuron is used for decoding freq&amp
%


%% perform population decoding

nIter = 200; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
minTrialNum = 15; %number of trials per condition to use
targetAmp = 5; % 0.33,0.76,1.6,3.5,[7.7]
targetFreq = 2; % 360, [525], 760, 1100, 1600
usePCA = true; % use PCA step prior to model training
kStart = 1; % first PC to include
kComp = 15; % total number of PCs to use in model training / testing
useFreqEffNeu = true; % use neurons with significant anova freq effects
useAllNeu = false; % use all available neurons for modeling
includeErrBar = true; %include error bars in plots
cellType = ["PV","SOM","EXC"]; % order for plotting
% cellType = "EXC"; % order for plotting
% selType = ["FreqMod","AmpMod","DualMod","InterMod"]; 
selType = "InterMod";
% is.allMod = (is.FreqMod | is.DualMod | is.InterMod);
% selType = "allMod";
% is.all = true(nCellTotal,1);
% selType = "all";
neuStepSize = 5; % num neurons to add to model each step
nCellDraws = 1; % over-estimation to poulate arrays
maxDraw = 100; % for EXC
matchedDraw = 25;
useCrossVal = false;

tic;

% CREATE OUTPUT ARRAYS
PCTall = zeros(2,1,nIter,length(selType),nCellDraws,length(cellType)); % (na,nj,nt,nd,nc)
MDLall = cell(2,1,nIter,length(selType),nCellDraws,length(cellType)); % (na,nj,nt,nd,nc)
Coefall = cell(2,1,nIter,length(selType),nCellDraws,length(cellType)); % (na,nj,nt,nd,nc)
CNall = cell(2,1,nIter,length(selType),nCellDraws,length(cellType)); % (na,nj,nt,nd,nc)

for nc = 1:length(cellType) %nt = 1 : length(selType)
disp(cellType(nc));
% disp(selType(nt));

for nt = 1 : length(selType) %nd = 1 : length(nCellDraws)
    fprintf(strcat(selType(nt)," "))
    %disp(strcat("Draw ",num2str(nCellDraws(nd))))

    % nCellDraws = kComp : neuStepSize : maxNeuAvail;

    % % MAX AVAIL. FOR PV AND SOM, 100 FOR EXC
    % if strcmp(cellType(nc),'EXC')
    %     nCellDraws = maxDraw;
    % else % PV | SOM
    %     nCellDraws = sum(is.(cellType(nc)) & is.(selType(nt)));
    % end
    % MATCHED DRAW SIZE FOR ALL CELL TYPES
    nCellDraws = matchedDraw;

    trainTrial = round(minTrialNum*trainPortion);
    testTrial = round(minTrialNum*(1-trainPortion));
    trainFreq = repmat(freqVal,trainTrial,1);
    testFreq = repmat(freqVal,testTrial,1);
    trainAmp = repmat(ampValAll,trainTrial,1);
    testAmp = repmat(ampValAll,testTrial,1);
    
    for nd = 1 : length(nCellDraws) %nc = 1:length(cellType)
        % fprintf(strcat(" ",cellType(nc)))
        % fprintf(strcat(num2str(nCellDraws(nd))," "));
        nCellDraw = nCellDraws(nd);

        % select cell type
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
            CNall{1,1,nj,nt,nd,nc} = selectedCells;
            CNall{2,1,nj,nt,nd,nc} = selectedCells;
           
            % pre-allocate training and testing data array
            trainRespFD = zeros(nFreq,trainTrial,nCellDraw);
            testRespFD = zeros(nFreq,testTrial,nCellDraw);
            trainRespAD = zeros(nFreq,trainTrial,nCellDraw);
            testRespAD = zeros(nFreq,testTrial,nCellDraw);
    
            % loop over selected cells
            for nk = 1:length(selectedCells)
                cellResp = Tc.responses{selectedCells(nk)};
                cellFreq = Tc.frequencies{selectedCells(nk)};
                cellAmp = Tc.amplitudes{selectedCells(nk)};
                          
                for nl = 1:nFreq
                    % randomly draw trials at each frequency, no replacement
                    trialNum = find(cellFreq == freqVal(nl) & cellAmp == ampValAll(targetAmp));
                    % trialNum = find(cellFreq == freqVal(nl)); % ALL AMPS
                    selectedTrials = randsample(trialNum, trainTrial+testTrial, false);
                    % populate training and testing arrays
                    trainRespFD(nl,:,nk) = cellResp(selectedTrials(1:trainTrial));
                    testRespFD(nl,:,nk) = cellResp(selectedTrials(trainTrial+1:end));
                end

                for nl = 1:nAmp
                    % randomly draw trials at each amplitude, no replacement
                    trialNum = find(cellAmp == ampValAll(nl) & cellFreq == freqVal(targetFreq));
                    % trialNum = find(cellAmp == ampValAll(nl)); % ALL FREQS
                    selectedTrials = randsample(trialNum, trainTrial+testTrial, false);
                    
                    % populate training and testing arrays
                    trainRespAD(nl,:,nk) = cellResp(selectedTrials(1:trainTrial));
                    testRespAD(nl,:,nk) = cellResp(selectedTrials(trainTrial+1:end));
                end
            end
    
            % Reshape matricies
            trainRespFD = reshape(trainRespFD, nFreq*trainTrial,nCellDraw);
            testRespFD = reshape(testRespFD, nFreq*testTrial,nCellDraw);
            trainRespAD = reshape(trainRespAD, nFreq*trainTrial,nCellDraw);
            testRespAD = reshape(testRespAD, nFreq*testTrial,nCellDraw);
    
            % Z-Scoring 
            for nz = 1 : nCellDraw
                mu = mean(trainRespFD(:,nz));
                sigma = std(trainRespFD(:,nz));
                trainRespFD(:,nz) = (trainRespFD(:,nz) - mu) ./ sigma;
                testRespFD(:,nz) = (testRespFD(:,nz) - mu) ./ sigma;

                mu = mean(trainRespAD(:,nz));
                sigma = std(trainRespAD(:,nz));
                trainRespAD(:,nz) = (trainRespAD(:,nz) - mu) ./ sigma;
                testRespAD(:,nz) = (testRespAD(:,nz) - mu) ./ sigma;
            end
    
            % FREQUENCY DECODING
            [coeff,score,latent] = pca(trainRespFD);
            explPerPC = latent ./ sum(latent) * 100;   % percent variance per PC
            cumExpl = cumsum(explPerPC); % cumulative percent explained
            trainResp_red = score(:,kStart:kComp);
            if useCrossVal == true % USE CROSS-VALIDATION
                mdlFD = fitcdiscr(trainResp_red,trainFreq,'DiscrimType','linear','CrossVal','on',...
                    'SaveMemory','off'); %train model, save all coeffs for each CV-fold
                PCTall(1,1,nj,nt,nd,nc) = 1 - kfoldLoss(mdlFD); %accuracy (1- classification error)
            else % NO CROSS-VALIDATION
                mdlFD = fitcdiscr(trainResp_red,trainFreq,'DiscrimType','linear'); %train model
                predFreq = predict(mdlFD,testRespFD*coeff(:,kStart:kComp));
                PCTall(1,1,nj,nt,nd,nc) = mean(predFreq == testFreq); % accuracy
            end
            % collect output
            MDLall{1,1,nj,nt,nd,nc} = mdlFD; % model parameters
            Coefall{1,1,nj,nt,nd,nc} = coeff(:,kStart:kComp);

            % AMPLITUDE DECODING
            [coeff,score,latent] = pca(trainRespAD);
            trainResp_red = score(:,kStart:kComp);
            if useCrossVal == true % USE CROSS-VALIDATION
                mdlAD = fitcdiscr(trainResp_red,trainAmp,'DiscrimType','linear','CrossVal','on',...
                    'SaveMemory','off'); %train model, save all coeffs for each CV-fold
                PCTall(2,1,nj,nt,nd,nc) = 1 - kfoldLoss(mdlAD); %accuracy (1- classification error)
            elseif useCrossVal == false % NO CROSS VALIDATION
                mdlAD = fitcdiscr(trainResp_red,trainAmp,'DiscrimType','linear'); %train model
                predAmp = predict(mdlAD,testRespAD*coeff(:,kStart:kComp));
                PCTall(2,1,nj,nt,nd,nc) = mean(predAmp == testAmp); % accuracy
            end
            % collect output
            MDLall{2,1,nj,nt,nd,nc} = mdlAD; % model parameters
            Coefall{2,1,nj,nt,nd,nc} = coeff(:,kStart:kComp);
        end 
    end
    fprintf(newline)
end

toc;

end





%% ACCURACY AT N=40, EACH OF 4 CELL TYPES

% targetAmp = 5; % 7.7 um
% xSpace = [1,2,3;6,7,8;11,12,13;16,17,18];
% ntOrder = [2 1 3 4];
% ntNames = {"A","F","F+A","F*A"};
% figure('theme','light','color',[1 1 1]); hold on;
% plot(-1,-1,'.','MarkerSize',8,'Color',colors.PV);
% plot(-1,-1,'.','MarkerSize',8,'Color',colors.SOM);
% plot(-1,-1,'.','MarkerSize',8,'Color',colors.EXC);
% for nt = 1 : length(selType)
%     for nc = 1 : length(cellType)
%         yVals = PCTall(2,1,:,ntOrder(nt),1,nc);
%         yMean = mean(yVals)*100; % mean y value fonverted to percentage
%         yCIs = (mkCI(yVals) - mean(yVals)).*100; % 95% confidence intervals centered at mean
%         plot(xSpace(nt,nc),yMean,'.','MarkerSize',8,'Color',colors.(cellType(nc)));
%         errorbar(xSpace(nt,nc),yMean,yCIs(1),yCIs(2),'Color',colors.(cellType(nc)));
%     end
% end
% 
% ylabel('model accuracy (%)','FontName','Arial');
% xlabel('Selectivity Type','FontName','Arial');
% title('Frequency decoding performance of each selectivity type','FontName','Arial');
% subtitle('Amplitude: 0.33\mum, model size: 40 neurons ea.','FontName','Arial');
% set(gca,'YLim',[15 65],'XLim',[0 18],'XTick',xSpace(:,2),'XTickLabel',ntNames,'FontName','Arial');
% legend({"PV","SOM","EXC"},'Location','Northwest','FontName','Arial','box','off');


%% RELATIONSHIP BETWEEN DELTA PREDICTOR [FD vs. AD]
% 
% targetCell = 1; % PV,SOM,EXC
% nr = 2; %targt iteration
% 
% mdlFD = MDLall{1,1,nr,1,1,targetCell};
% mdlAD = MDLall{2,1,nr,1,1,targetCell};
% figure('theme','light','color',[1 1 1]); hold on;
% plot(mdlFD.DeltaPredictor,mdlAD.DeltaPredictor,'.','Color',colors.PV);

%% VISUALIZE decoder weights

% for nt = 1 : 2
%     figure('theme','light','color',[1 1 1]); hold on;
%     % custom color map for weights display
%     vals = [-6, 0, 6];
%     cols = [0 0 1; 0.5 0.5 0.5; 1 0 0];
%     xi = linspace(-6,6,256);
%     cmap = interp1(vals, cols, xi, 'linear');
%     colormap(cmap);
%     mdl = MDLall{nt,1,1,1,1,1};
%     pos = [2 3 4 5 8 9 10 14 15 20];
%     indx = 1;
%     for n1 = 1 : 5
%         for n2 = n1+1 : 5
%             locMat = reshape(mdl.Coeffs(n1,n2).Linear,5,3);
%             % locMat = reshape(mdlAD.Coeffs(n1,n2).Linear,5,3);
%             subplot(5,5,pos(indx));
%             imagesc(locMat);
%             if nt == 1
%                 title(strcat(num2str(freqVal(n1))," v. ",num2str(freqVal(n2))));
%             elseif nt == 2
%                 title(strcat(num2str(ampValAll(n1))," v. ",num2str(ampValAll(n2))));
%             end
%             clim([-6 6]);
%             set(gca,'XColor',[1 1 1],'YColor',[1,1,1])
%             indx = indx + 1;
%         end
%     end
% end

%% COMPARISON OF ACCURACY, WITH & WITHOUT CROSSVAL
% cd('C:\Users\skich\Desktop\WORK');
% load('modelResults_ampAndFreq_decodeSequentially_CrossV_ExternalV.mat');

% PCTall_CV = squeeze(PCTall_CV);
% PCTall_EV = squeeze(PCTall_EV);
% 
% xSpace = [1,2,3; 5,6,7];
% 
% nd = 2; % decoding type [ FD | AD ]
% nt = 4; % selType [ F | A | F+A | F*A ]
% 
% figure('theme','light','color',[1 1 1]); hold on;
% plot(-1,-1,'.','Color',colors.PV);
% plot(-1,-1,'.','Color',colors.SOM);
% plot(-1,-1,'.','Color',colors.EXC);
% 
% for nc = 1 : 3
%     % EXTERNAL VALIDATION
%     yMean = mean(PCTall_EV(nd,:,nt,nc));
%     yCI = mkCI(squeeze(PCTall_EV(nd,:,nt,nc)));
%     errorbar(xSpace(1,nc),yMean,yMean-yCI(1),yMean-yCI(2),'Color',colors.(cellType(nc)));
%     plot(xSpace(1,nc),yMean,'.','Color',colors.(cellType(nc)));
% 
%     % CROSS VALIDATION
%     yMean = mean(PCTall_CV(nd,:,nt,nc));
%     yCI = mkCI(squeeze(PCTall_CV(nd,:,nt,nc)));
%     errorbar(xSpace(2,nc),yMean,yMean-yCI(1),yMean-yCI(2),'Color',colors.(cellType(nc)));
%     plot(xSpace(2,nc),yMean,'.','Color',colors.(cellType(nc)));
% end
% set(gca,'FontName','Arial','XTick',xSpace(:,2),'XTickLabels',{"Ext-Validation","Cross-Validation"},...
%     'XLim',[0 max(xSpace,[],"all")+1]);
% if nd == 1 
%     title('Freq-decode performance by validation type','FontName','Arial');
%     ylabel('Frequency decoding accuracy %','FontName','Arial');
% elseif nd == 2
%     title('Amp-decode performance by validation type','FontName','Arial');
%     ylabel('Amplitude decoding accuracy %','FontName','Arial');
% end
% legend({"PV","SOM","EXC"},'Box','off','location','northeast')


%% COMPARISON OF ACCURACY, WITH & WITHOUT CROSSVAL
% ALL AMPLITUDES TESTED
% cd('C:\Users\skich\Desktop\WORK');
% load('modelResults_ampAndFreq_decodeSequentially_5x5.mat');

% PCTall5x5 = squeeze(PCTall5x5);
% 
% xSpace = [1,2,3; 5,6,7];
% 
% nd = 1; % decoding type [ FD | AD ]
% nt = 4; % selType [ F | A | F+A | F*A ]
% na = 5; % amplitude
% 
% figure('theme','light','color',[1 1 1]); hold on;
% plot(-1,-1,'.','Color',colors.PV);
% plot(-1,-1,'.','Color',colors.SOM);
% plot(-1,-1,'.','Color',colors.EXC);
% 
% for nc = 1 : 3
%     % EXTERNAL VALIDATION
%     yMean = mean(PCTall5x5(nd,:,nt,nc,na));
%     yCI = mkCI(squeeze(PCTall5x5(nd,:,nt,nc,na)));
%     errorbar(xSpace(1,nc),yMean,yMean-yCI(1),yMean-yCI(2),'Color',colors.(cellType(nc)));
%     plot(xSpace(1,nc),yMean,'.','Color',colors.(cellType(nc)));
% 
%     % CROSS VALIDATION
%     yMean = mean(PCTall5x5(nd,:,nt,nc));
%     yCI = mkCI(squeeze(PCTall5x5(nd,:,nt,nc)));
%     errorbar(xSpace(2,nc),yMean,yMean-yCI(1),yMean-yCI(2),'Color',colors.(cellType(nc)));
%     plot(xSpace(2,nc),yMean,'.','Color',colors.(cellType(nc)));
% end
% set(gca,'FontName','Arial','XTick',xSpace(:,2),'XTickLabels',{"Ext-Validation","Cross-Validation"},...
%     'XLim',[0 max(xSpace,[],"all")+1]);
% if nd == 1 
%     title('Freq-decode performance by validation type','FontName','Arial');
%     subtitle(strcat("Amplitude: ",num2str(ampValAll(na)),"\mum"),'FontName','Arial');
%     ylabel('Frequency decoding accuracy %','FontName','Arial');
% elseif nd == 2
%     title('Amp-decode performance by validation type','FontName','Arial');
%     subtitle(strcat("Frequency: ",num2str(freqVal(na)),"Hz"),'FontName','Arial');
%     ylabel('Amplitude decoding accuracy %','FontName','Arial');
% end
% legend({"PV","SOM","EXC"},'Box','off','location','northeast')

%% targeted analysis of neural correlates

% cd('C:\Users\skich\Desktop\WORK');
% load('modelResults_ampAndFreq_decodeSequentially.mat');

targetSel = 1; % F, A, F+A, F*A
% targetRhoVals = [0.41, 0.48, 0.49];
matRP = zeros(length(cellType),nIter,2);
weights = struct();
for nc = 1 : length(cellType)
    % output matrix
    coefFD = Coefall{1,1,1,targetSel,1,nc};
    neuMat = zeros(nIter,size(coefFD,1),10,2); % [Iterations] x [neurons] x [boundary] x [FD|AD]
    for nr = 1 : nIter
        % extract variables from rep
        mdlFD = MDLall{1,1,nr,targetSel,1,nc};
        coefFD = Coefall{1,1,nr,targetSel,1,nc};
        mdlAD = MDLall{2,1,nr,targetSel,1,nc};
        coefAD = Coefall{2,1,nr,targetSel,1,nc};
        % extract FD weights
        indx = 1;
        for n1 = 1 : 5
            for n2 = n1+1 : 5
                PCWeights = mdlFD.Coeffs(n1,n2).Linear; % linear weights for discriminating n1 v. n2
                neuWeights = coefFD * PCWeights; % PC loading matrix (40x15) * LDA weights PC space (15x1)
                neuMat(nr,:,indx,1) = neuWeights;
                indx = indx + 1;
            end
        end 
        % extract AD weights
        indx = 1;
        for n1 = 1 : 5
            for n2 = n1+1 : 5
                PCWeights = mdlAD.Coeffs(n1,n2).Linear; % linear weights for discriminating n1 v. n2
                neuWeights = coefAD * PCWeights; % PC loading matrix (40x15) * LDA weights PC space (15x1)
                neuMat(nr,:,indx,2) = neuWeights;
                indx = indx + 1;
            end
        end
        FDweights = max(abs(squeeze(neuMat(nr,:,:,1))),[],2);
        ADweights = max(abs(squeeze(neuMat(nr,:,:,2))),[],2);
        [rho,pval] = corr(cat(2,FDweights,ADweights),'Type','Spearman');
        matRP(nc,nr,:) = [rho(1,2) pval(1,2)];
    end
    weights.(cellType(nc)) = neuMat; % nIter(200) x nNeurons(25) x nBoundaries(10) x FD-or-AD(2)
end

% % PLOT AVG RHO BY CELL TYPE
% figure('theme','light','color',[1 1 1]); hold on;
% for nc = 1 : 3
%     % locRhoVals = matRP(nc,matRP(nc,:,2)<0.05,1); % significant rho values only
%     % locRhoVals = locRhoVals(1:200); %isolate 200 sig reps
%     locRhoVals = matRP(nc,:,1); % all rho values
%     locCIvals = mkCI(locRhoVals) - mean(locRhoVals);
%     errorbar(nc,mean(locRhoVals),locCIvals(1),locCIvals(2),'Color',colors.(cellType(nc)));
%     plot(nc,mean(locRhoVals),'.','MarkerSize',10,'Color',colors.(cellType(nc)));
% end
% ylabel("Spearman's rho",'FontName','Arial');
% title("Correlation between neuron decoding weights",'FontName','Arial');
% subtitle("Frequency vs. Amplitude models (n=200)",'FontName','Arial');
% set(gca,'XLim',[0.5 3.5],'XTick',1:3,'XTickLabel',["PV","SOM","EXC"],...
%     'FontName','Arial');
    

% PLOT SCATTER MAX WEIGHT PER NEURON FROM EXAMPLE MODELS
for nc = 1 : 3
    % EXTRACT VALUES
    % locRhoVals = matRP(nc,matRP(nc,:,2)<0.05,1); % significant rho values only
    locRhoVals = matRP(nc,:,1); % all rho values
    % locRhoVals = locRhoVals(1:200); %isolate 200 sig reps

    targetRho = mean(locRhoVals); 
    % targetRho = mean(matRP(nc,:,1)); % all rho values
    [~,indx] = sort(abs(matRP(nc,:,1)-targetRho),'ascend'); % find model rep closest to average rho
    %[~,indx] = sort(abs(matRP(nc,:,1)-targetRhoVals(nc)),'ascend'); % find model rep closest specified rho (above)
    % [~,indx] = max(matRP(nc,:,1)); % use model rep with highest rho
    neuMat = weights.(cellType(nc));
    FDweights = max(abs(squeeze(neuMat(indx(1),:,:,1))),[],2);
    ADweights = max(abs(squeeze(neuMat(indx(1),:,:,2))),[],2);

    % PLOTTING
    figure('theme','light','color',[1 1 1]); hold on;
    plot(FDweights,ADweights,'.','Color',colors.(cellType(nc)),'MarkerSize',8);
    xlabel("Max neuron weight (frequency decode)",'FontName','Arial');
    ylabel("Max neuron weight (amplitude decode)",'FontName','Arial');
    titleTxt = strcat(cellType(nc)," Neuron use in each decoding condition");
    title(titleTxt,'FontName','Arial');
    % stitleTxt = strcat("rho = ",num2str(matRP(nc,indx(1),1),3),"   pval = ",num2str(matRP(nc,indx(1),2),1)); % rho and p-val
    stitleTxt = strcat("Spearman's rho (\rho) = ",num2str(matRP(nc,indx(1),1),3)); % rho and p-val
    subtitle(stitleTxt,'FontName','Arial');
    locXLim = [min(FDweights)-.1 max(FDweights)+.1];
    locYLim = [min(ADweights)-.1 max(ADweights)+.1];
    set(gca,'XLim',locXLim,'XLimMode','manual','YLim',locYLim,'YLimMode','manual');

    % LINEAR FIT
    lmdl = fit(FDweights,ADweights,'poly1');
    plot(locXLim,lmdl(locXLim),'--','Color',colors.(cellType(nc)));

    % % FIT TEXT
    % s0 = "Linear fit: - - -";
    % s1 = strcat("f(x) = ",num2str(lmdl.p1,3),"x + ",num2str(lmdl.p2,3));
    % s2 = strcat("p1 bounds: [",num2str(ci(1,1),3)," ",num2str(ci(2,1),3),"]");
    % s3 = strcat("p2 bounds: [",num2str(ci(1,2),3)," ",num2str(ci(2,2),3),"]");
    % text(3, 0.4, {s0, s1, s2, s3},'HorizontalAlignment','right','VerticalAlignment','middle',...
    %     'FontSize',10,'FontName','Arial','Color',colors.(cellType(nc)));
    % hold off;
end

%% scatter plot of significant rho values by each cell type
figure('theme','light','color',[1 1 1]); hold on;
colorMatPV = double(matRP(1,:,2) < 0.05)' + 1;
colorMatSOM = double(matRP(2,:,2) < 0.05)' + 3;
colorMatEXC = double(matRP(3,:,2) < 0.05)' + 5;
colorMat = cat(1,colorMatPV, colorMatSOM, colorMatEXC);
allXVals = cat(1,ones(nIter,1),ones(nIter,1).*2,ones(nIter,1).*3);
allYVals = cat(1,squeeze(matRP(1,:,1))',squeeze(matRP(2,:,1))',squeeze(matRP(3,:,1))');
tbl = table(allXVals,allYVals,colorMat,'VariableNames',{'X','Y','Color'});
swarmchart(tbl,'X','Y','filled','ColorVariable','Color');
ylabel("Spearman's Rho",'FontName','Arial');
xlabel("",'FontName','Arial');
title("Rho values per cell type over 200 repetitions",'FontName','Arial');
set(gca,'XTick',1:3,'XTickLabel',cellType,'FontName','Arial');
map = [colors.PV+.1; colors.PV-.1; colors.SOM+.1; colors.SOM-.1; colors.EXC+.1; colors.EXC-.1];
colormap(map)

%% COMPARE RHO VALUES: FISHER Z-TRANSFORMATION & REPEATED MEASURES ANOVA


% rho = zeros(3,200);
% for nc = 1 : 3
%     % EXTRACT VALUES
%     locRhoVals = matRP(nc,matRP(nc,:,2)<0.05,1); % significant rho values only
%     rho(nc,:) = locRhoVals(1:200); %isolate 200 sig reps


% % REPEATED MEASURES ANOVA
% % Build table and within-design for repated measures anova
% T = array2table(z.','VariableNames',{'Cond1','Cond2','Cond3'});
% within = table(categorical({'C1';'C2';'C3'}),'VariableNames', {'Condition'});
% % Fit model
% rm = fitrm(T, 'Cond1-Cond3 ~ 1', 'WithinDesign', within);
% ranovatbl = ranova(rm);
% disp(ranovatbl);
% % Pairwise Comparisons
% mc = multcompare(rm, 'Condition', 'ComparisonType', 'bonferroni');
% disp(mc);

% PERMUTATION TEST: null hypothesis = rho centered at zero
nPerm = 10000; % permutations
for nc = 1 : 3 % cellType
    rho = matRP(:,:,1); 
    z = atanh(rho);
    locRho = z(nc,:)'; % 200×1 vector
    % Sign-flip permutations
    permStats = zeros(nPerm,1);
    for i = 1:nPerm
        signs = sign(rand(size(locRho)) - 0.5);
        permStats(i) = mean(locRho .* signs);
    end
    % One-sided test: above zero
    pval = mean(permStats >= mean(locRho));
end


% PLOT SWARM CHART FISHER Z-TRANSFORM SPEARMAN RHO
figure('theme','light','color',[1 1 1]); hold on;
for nc = 1 : 3
    swarmchart(ones(1,nIter).*nc,z(nc,:),8,'filled','o','XJitterWidth',0.5,'MarkerFaceColor',colors.(cellType(nc)),...
        'MarkerFaceAlpha',0.9);
    % % 95 % CIs
    % locCI = mkCI(z(nc,:));
    % plot([-0.5 0.5]+nc,[locCI(1) locCI(1)],'--','Color',colors.(cellType(nc)));
    % plot([-0.5 0.5]+nc,[locCI(2) locCI(2)],'--','Color',colors.(cellType(nc)));
    % 2*Sigma bounds
    locSig = std(z(nc,:));
    locBounds = median(z(nc,:)) + (2*locSig*[-1 1]);
    plot([-0.5 0.5]+nc,repmat(locBounds(1),1,2),'--','Color',colors.(cellType(nc)));
    plot([-0.5 0.5]+nc,repmat(locBounds(2),1,2),'--','Color',colors.(cellType(nc)));
end
set(gca,'XLim',[0.5 3.5],'FontName','Arial','XTick',1:3,'XTickLabel',cellType);
ylabel("Fisher z-transformed Spearman's \rho",'FontName','Arial');
title("Correlations of neuron use in each decoding condition",'FontName','Arial');
subtitle("200 repetitions per type, freq-decode vs. amp-decode",'FontName','Arial');

% % PLOT STATS MULTIPLE COMPARISONS TEST
% pValThresh = 0.05;
% maxY = repmat(max(max(z))+0.02,1,2);
% if mc.pValue(1) < pValThresh % PV vs SOM
%     plot([1.1 1.9],maxY,'k-','LineWidth',1.2);
%     text(1.5,maxY(1),"**",'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
% end
% if mc.pValue(2) < pValThresh % PV vs EXC
%     plot([1.1 2.9],maxY+0.04,'k-','LineWidth',1.2);
%     text(2,maxY(1)+0.04,"**",'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
% end
% if mc.pValue(4) < pValThresh % SOM vs EXC
%     plot([2.1 2.9],maxY,'k-','LineWidth',1.2);
%     text(2.5,maxY(1),"*",'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
% end

% PERMUTATION TESTS: null hypothesis = no pairwise difference between groups
maxY = repmat(max(max(z))+0.02,1,2);
idxOrd = [1 2; 1 3; 2 3]; % PV v SOM; PV v EXC; SOM v EXC
pval = zeros(size(idxOrd,1),1);
for nc = 1 : size(idxOrd,1)
    x = z(idxOrd(nc,1),:)';
    y = z(idxOrd(nc,2),:)';
    combined = [x; y];
    permStats = zeros(nPerm,1);
    for i = 1:nPerm
        idx = randperm(numel(combined));
        permX = combined(idx(1:numel(x)));
        permY = combined(idx(numel(x)+1:end));
        permStats(i) = mean(permX) - mean(permY);
    end
    pval(nc) = mean(abs(permStats) >= abs(mean(x) - mean(y)));
end

% PLOT significance
maxY = repmat(max(max(z))+0.02,1,2);
if pval(1) < pValThresh % PV vs SOM
    plot([1.1 1.9],maxY,'k-','LineWidth',1.2);
    text(1.5,maxY(1),"**",'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
end
if pval(2) < pValThresh % PV vs EXC
    plot([1.1 2.9],maxY+0.04,'k-','LineWidth',1.2);
    text(2,maxY(1)+0.04,"**",'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
end
if pval(3) < pValThresh % SOM vs EXC
    % plot([2.1 2.9],maxY,'k-','LineWidth',1.2);
    % text(2.5,maxY(1),"*",'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
end



%% 
% % correlation between avg FD vs AD weight from each neuron
% [rho,pval] = corr(cat(2, mean(abs(neuMat(:,:,1)),2), mean(abs(neuMat(:,:,2)),2)),'Type','Spearman');

% correlation between max FD weight vs AD weight from each neuron
locMat = cat(2,max(abs(neuMat(:,:,1)),[],2),max(abs(neuMat(:,:,2)),[],2));
[rho,pval] = corr(locMat,'Type','Spearman');

figure('theme','light','color',[1 1 1]); hold on;
plot(locMat(:,1),locMat(:,2),'.','Color',colors.(cellType(targetCell)),'MarkerSize',8);
set(gca,'XLim',[0 3.5],'YLim',[0 3]);
xlabel("Max neuron weight (frequency decode)",'FontName','Arial');
ylabel("Max neuron weight (amplitude decode)",'FontName','Arial');
titleTxt = strcat(cellType(targetCell)," Neuron use in each decoding condition");
title(titleTxt,'FontName','Arial');
stitleTxt = strcat("rho = ",num2str(rho(1,2),3),"   pval = ",num2str(pval(1,2),1));
subtitle(stitleTxt,'FontName','Arial');

% Linear fit
lmdl = fit(locMat(:,1),locMat(:,2),'poly1');
ci = confint(lmdl); 
plot([0 3.5],lmdl([0 3.5]),'--','Color',colors.(cellType(targetCell)));
% fit text
s1 = strcat("f(x) = ",num2str(lmdl.p1,3),"x + ",num2str(lmdl.p2,3));
s2 = strcat("p1 bounds: [",num2str(ci(1,1),3)," ",num2str(ci(2,1),3),"]");
s3 = strcat("p2 bounds: [",num2str(ci(1,2),3)," ",num2str(ci(2,2),3),"]");
text(3.5, 0.4, {s1, s2, s3},'HorizontalAlignment','right','VerticalAlignment','middle',...
    'FontSize',10,'FontName','Arial','Color',colors.(cellType(targetCell)));
hold off;



%% CORRELATION (SPEARMAN RHO) OF ALL SIGNIFICANT REPS ACROSS CELL TYPE and ONE CONDITON
targetSel = 1; % only F*A
% neuMatFD = nan(height(Tc),nIter);
matRP = zeros(200,2,3);
for nc = 1 : 3
    targetCell = nc; % PV,SOM,EXC
    for nr = 1 : nIter % loop model iterations
        % extract variables from rep
        mdlFD = MDLall{1,1,nr,targetSel,1,targetCell};
        coefFD = Coefall{1,1,nr,targetSel,1,targetCell};
        mdlAD = MDLall{2,1,nr,targetSel,1,targetCell};
        coefAD = Coefall{2,1,nr,targetSel,1,targetCell};
        cellNums = CNall{1,1,nr,targetSel,1,targetCell};
        
        % output matrix
        neuMat = zeros(size(coefAD,1),10,2); % [neurons] x [boundary] x [FD|AD]
        
        % extract FD weights
        indx = 1;
        for n1 = 1 : 5
            for n2 = n1+1 : 5
                PCWeights = mdlFD.Coeffs(n1,n2).Linear; % linear weights for discriminating n1 v. n2
                neuWeights = coefFD * PCWeights; % PC loading matrix (40x15) * LDA weights PC space (15x1)
                neuMat(:,indx,1) = neuWeights;
                indx = indx + 1;
            end
        end 
    
        % extract AD weights
        indx = 1;
        for n1 = 1 : 5
            for n2 = n1+1 : 5
                PCWeights = mdlAD.Coeffs(n1,n2).Linear; % linear weights for discriminating n1 v. n2
                neuWeights = coefAD * PCWeights; % PC loading matrix (40x15) * LDA weights PC space (15x1)
                neuMat(:,indx,2) = neuWeights;
                indx = indx + 1;
            end
        end 
    
        % correlation between avg FD vs AD weight from each neuron
        locMat = cat(2, max(abs(neuMat(:,:,1)),[],2), max(abs(neuMat(:,:,2)),[],2));
        [rho,pval] = corr(locMat,'Type','Spearman');
        matRP(nr,1,nc) = rho(1,2);
        matRP(nr,2,nc) = pval(1,2);
    
    end
end


figure('theme','light','color',[1 1 1]); hold on;

for nc = 1 : 3
    locRhoVals = matRP(matRP(:,2,nc) < 0.05,1,nc);
    locCIvals = mkCI(locRhoVals) - mean(locRhoVals);
    errorbar(nc,mean(locRhoVals),locCIvals(1),locCIvals(2),'Color',colors.(cellType(nc)));
    plot(nc,mean(locRhoVals),'.','MarkerSize',10,'Color',colors.(cellType(nc)));
end
ylabel("Spearman's rho",'FontName','Arial');
title("Correlation between neuron decoding weights",'FontName','Arial');
subtitle("Frequency vs. Amplitude models (n=200)",'FontName','Arial');
set(gca,'YLim',[0.26 0.44],'XLim',[0.5 3.5],'XTick',1:3,'XTickLabel',["PV","SOM","EXC"],...
    'FontName','Arial');



%% targeted analysis of individual rep
targetSel = 4; % F, A, F+A, F*A
targetCell = 3; % PV,SOM,EXC
% neuMatFD = nan(height(Tc),nIter);
matRP = zeros(200,2);
nr = 1;

%for nr = 1 : nIter % loop model iterations

% extract variables from rep
mdlFD = MDLall{1,1,nr,targetSel,1,targetCell};
coefFD = Coefall{1,1,nr,targetSel,1,targetCell};
mdlAD = MDLall{2,1,nr,targetSel,1,targetCell};
coefAD = Coefall{2,1,nr,targetSel,1,targetCell};
cellNums = CNall{1,1,nr,targetSel,1,targetCell};

% output matrix
neuMat = zeros(40,10,2); % [neurons] x [boundary] x [FD|AD]

% extract FD weights
indx = 1;
for n1 = 1 : 5
    for n2 = n1+1 : 5
        PCWeights = mdlFD.Coeffs(n1,n2).Linear; % linear weights for discriminating n1 v. n2
        neuWeights = coefFD * PCWeights; % PC loading matrix (40x15) * LDA weights PC space (15x1)
        neuMat(:,indx,1) = neuWeights;
        indx = indx + 1;
    end
end 

% extract AD weights
indx = 1;
for n1 = 1 : 5
    for n2 = n1+1 : 5
        PCWeights = mdlAD.Coeffs(n1,n2).Linear; % linear weights for discriminating n1 v. n2
        neuWeights = coefAD * PCWeights; % PC loading matrix (40x15) * LDA weights PC space (15x1)
        neuMat(:,indx,2) = neuWeights;
        indx = indx + 1;
    end
end 

% correlation between avg FD vs AD weight from each neuron
[rho,pval] = corr(cat(2, mean(abs(neuMat(:,:,1)),2), mean(abs(neuMat(:,:,2)),2)),'Type','Spearman');
matRP(nr,1) = rho(1,2);
matRP(nr,2) = pval(1,2);

figure('theme','light','color',[1 1 1]); hold on;
plot(mean(abs(neuMat(:,:,1)),2),mean(abs(neuMat(:,:,2)),2),'.','Color',colors.(cellType(targetCell)));
xlabel("Avg. weight frequency decoding",'FontName','Arial');
ylabel("Avg. weight amplitude decoding",'FontName','Arial');
titleTxt = strcat(cellType(targetCell)," Neuron use in each decoding condition");
title(titleTxt,'FontName','Arial');
stitleTxt = strcat("rsq = ",num2str(rho(1,2)^2),"   pval = ",num2str(pval(1,2)));
subtitle(stitleTxt,'FontName','Arial');

%end


figure('theme','light','color',[1 1 1]); hold on;
% custom color map for weights display
vals = [-4, 0, 4];
cols = [0 0 1; 0.5 0.5 0.5; 1 0 0];
xi = linspace(-4,4,256);
cmap = interp1(vals, cols, xi, 'linear');
colormap(cmap);
tiledlayout(1,2);
nexttile(1)
imagesc(neuMat(:,:,1));
nexttile(2);
imagesc(neuMat(:,:,2));















