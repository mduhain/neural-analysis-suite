% FAS_multiCompare
%
%
%

cd('C:\Users\skich\Desktop\WORK');
% load("modelResults_dataset1_multiCompare.mat");

figure('theme','light','color',[1 1 1]); hold on;


% BASELINE
xSp = 1:3;
for nc = 1 : 3
    yMean = mean(PCTall_baseline(:,nc)).*100;
    yCIs = mkCI(PCTall_baseline(:,nc)).*100;
    plot([xSp(nc) xSp(nc)],yCIs,'-','LineWidth',1.5,'Color',colors.(cellType(nc)));
    plot(xSp(nc),yMean,'.','MarkerSize',12,'Color',colors.(cellType(nc)));
end

% % NO PC1
% xSp = 5:7;
% for nc = 1 : 3
%     yMean = mean(PCTall_noPC1(:,nc));
%     yCIs = mkCI(PCTall_noPC1(:,nc));
%     plot([xSp(nc) xSp(nc)],yCIs,'-','LineWidth',1,'Color',colors.(cellType(nc)));
%     plot(xSp(nc),yMean,'.','MarkerSize',8,'Color',colors.(cellType(nc)));
% end

% no rsq [0.9-1.0]
xSp = 6:8;
for nc = 1 : 3
    yMean = mean(PCTall_noAmpCorr(:,nc,2)).*100;
    yCIs = mkCI(PCTall_noAmpCorr(:,nc,2)).*100;
    plot([xSp(nc) xSp(nc)],yCIs,'-','LineWidth',1.5,'Color',colors.(cellType(nc)));
    plot(xSp(nc),yMean,'.','MarkerSize',12,'Color',colors.(cellType(nc)));
end

% no rsq [0.8-1.0]
xSp = 11:13;
for nc = 1 : 3
    yMean = mean(PCTall_noAmpCorr(:,nc,3)).*100;
    yCIs = mkCI(PCTall_noAmpCorr(:,nc,3)).*100;
    plot([xSp(nc) xSp(nc)],yCIs,'-','LineWidth',1.5,'Color',colors.(cellType(nc)));
    plot(xSp(nc),yMean,'.','MarkerSize',12,'Color',colors.(cellType(nc)));
end

hold off;

set(gca,'XLim',[0 14],'XTick',[2 7 12],'XTickLabel',{"Baseline","<0.9 amp-rsq.","<0.8 amp-rsq."},...
    'XTickLabelRotation',0,'YLim',[77 88],'YTick',78:2:88,'FontName','Arial');
ylabel('model accuracy (%)','FontName','Arial');
title('Frequency decoding performance by cell type','FontName','Arial');
subtitle('Remove neurons with responses correlated to amplitude','FontName','Arial');
legend({"PV","","SOM","","EXC"},'Location','NorthEast');