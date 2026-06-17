% tuningCurveAF.m
%   [] = tuningCurveAF(respStruct)
%   Accepts a response struct from R{} with amp/freq modulation
%   Plots amp-colored tuning curve in current figure

function [] = tuningCurveAF(respStruct)

    freqs = [360 525 760 1100 1600];
    mat = respStruct.avgRespPost;
    locRespVect = respStruct.avgRespAllTrials(:);
    indxSpaceMat = imrotate([1:5;21:25;41:45;61:65;81:85],90);
    indxSpace = indxSpaceMat(:); %vectorize
    specColor = [1,0,0; 0.6,0,0; 0.4,0.4,0.4; 0,0.6,0.6; 0,1,1]; % RED TO CYAN
    % specColor = [0,0,0; 0.2 0.2 0.2; 0.4 0.4 0.4; 0.6,0.6,0.6; 0.8,0.8,0.8]; % GRAYSCALE
    locColors = repmat(specColor,5,1);
    hold on;
    for ns = 1 : length(indxSpace)
        ymax = mean(locRespVect{ns}) + sem(locRespVect{ns});
        ymin = mean(locRespVect{ns}) - sem(locRespVect{ns});
        plot([indxSpace(ns) indxSpace(ns)],[ymax ymin],'LineStyle','-','Color',...
            locColors(ns,:),'LineWidth',1.6);
        plot(indxSpace(ns),mean(locRespVect{ns}),'.','Color',...
            locColors(ns,:),'MarkerSize',10);
    end
    for n = 1 : 5
        plot(indxSpaceMat(n,:),mat(n,:),'LineStyle','-','LineWidth',1.6,...
            'Color',locColors(n,:));
    end
    % adjust figure properties
    set(gca,'XTick',indxSpaceMat(3,:))
    set(gca,'XTickLabel',freqs);
    set(gca,'XLim',[min(indxSpace)-4 max(indxSpace)+4]);
    %ylabel('dF/F'); xlabel('Frequency (Hz)');
    %title("Frequency Tuning by Amplitude");
    %subtitle(strcat("Example Neuron : ",num2str(neuNum)));
end



