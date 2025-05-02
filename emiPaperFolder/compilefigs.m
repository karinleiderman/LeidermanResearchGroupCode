% clear all
fig1 = openfig('KDX_KDIXa.fig', 'reuse');
fig2 = openfig('KDX_kcat.fig', 'reuse');
fig3 = openfig('KDIXa_kcat.fig', 'reuse');

%load('AlphaBetaFitting/custom_cmap.mat')

figure
clf
t = tiledlayout(2, 2); % 1 row, 3 columns

ax1 = gca(fig1);
newAx1 = nexttile(t);
text(0.025,0.95,chars{1},'Units','normalized','Color',[1,1,1],'FontSize',labelFont,'FontWeight','bold')
copyAxesProperties(ax1, newAx1);
%colormap(custom_cmap)
Kd_IXaM=5.5276e+03;
Kd_MX=55.7235;
marker = 58;
hold on
red=[1 0 0];
color1 = red; % Pure red
color2 = red + (1 - red) * 0.25; % Red with 33% white
color3 = red + (1 - red) * 0.5; % Red with 66% white
color4 = red + (1 - red) * 0.75; % Red with 66% white
%color4 = [1, 1, 1]; % Pure white

% color1=[0 0 0];
% color2=[0.5 0.5 0.5];
% color3=[0.7 0.7 0.7];
% color4=[0.9 0.9 0.9];
plot(Kd_IXaM,Kd_MX,'.','color',color1,MarkerSize=marker)
plot(Kd_IXaM/10,Kd_MX,'.','color',color2,MarkerSize=marker)
plot(Kd_IXaM/100,Kd_MX,'.','color',color3,MarkerSize=marker)
plot(Kd_IXaM/1000,Kd_MX,'.','color',color4,MarkerSize=marker)
caxis([0 0.125])
xlabel('K_D^{IXa,bsAb} (nM)')
ylabel('K_D^{X,bsAb} (nM)')
set(gca,'FontSize',20)

yticks([1 100 10000]);
yticklabels({'10^0','10^2','10^4'})
xticks([1 100 10000]);
xticklabels({'10^0','10^2','10^4'})

%%%%%%%%%%%%%%%%%%%%%

newAx2=nexttile(t);


%DB(1,:) = load('DSB_IXa_pt0001Emi.dat');
DB(1,:) = load('DSB_IXa_pt001Emi.dat');
DB(2,:) = load('DSB_IXa_pt01Emi.dat');
DB(3,:) = load('DSB_IXa_pt1Emi.dat');
DB(4,:) = load('DSB_IXa_Emi.dat');
DB(5,:) = load('DSB_IXa_10Emi.dat');
DB(6,:) = load('DSB_IXa_100Emi.dat');


% 1 IXa in solution         IXa(pp,i)./IXatot,...
% 2 IXa on surface          IXa_bp(pp,i)./IXatot,...
% 3 IXa:X on surface        CZ(pp,i)./IXatot,...
% 4 IXa:M in solution       IXaM(pp,i)./IXatot,...
% 5 IXa:M on surface        IXaM_bp(pp,i)./IXatot,...
% 6 IXa:M:X on surface      IXaMX_bp(pp,i)./IXatot,...
% 7 IXa:M:Xa on surface     IXaMXa_bp(pp,i)./IXatot];


order = [1,2,3,7,6,5,4];
% lorder =['IXa in solution'...
%     'IXa on surface',...
%     'IXa:FX on surface',...
%      'IXa:M:FXa on surface',...
%      'IXa:M:FX on surface',...
%      'IXa:M on surface',...
%      'IXa:M in solution']


for i=1:6
    %data(i,:) = [sum(DB(i,[1:3,7])) DB(i,[6,5,4])];
    data(i,:)=DB(i,order);
end


 % X = {'0.01\cdotK_{D,emi}^{IXa,M}',...
 %     '0.1\cdotK_{D,emi}^{IXa,M}','K_{D,emi}^{IXa,M}',...
 %     '10\cdotK_{D,emi}^{IXa,M}','100\cdotK_{D,emi}^{IXa,M}'};
 % 
 % X = {'0.01\timesK_{D,emi}^{IXa,M}',...
 %     '0.1\timesK_{D,emi}^{IXa,M}','K_{D,emi}^{IXa,M}',...
 %     '10\timesK_{D,emi}^{IXa,M}','100\timesK_{D,emi}^{IXa,M}'};

  X = {'0.001\timesK_{D}^{IXa,M}',...
      '0.01\timesK_{D}^{IXa,M}',...
     '0.1\timesK_{D}^{IXa,M}','K_{D}^{IXa,M}'};

b=bar(X,data(1:4,:),'stacked')
set(b,'Facecolor','Flat')

cmap=colormap(newAx2,parula(7));
%cmap2=colormap(newAx3,parula(7));
% for i=1:7
%     b(i).CData = cmap(i,:);
% end
% 
b(1).CData = cmap(3,:);
b(2).CData = cmap(4,:);
b(3).CData = cmap(5,:);
b(4).CData = cmap(6,:);
b(5).CData = cmap(7,:);
b(6).CData = cmap(1,:);
b(7).CData = cmap(2,:);

% 
% b(1).CData = cmap(4,:);
% b(2).CData = cmap(3,:);
% b(3).CData = cmap(2,:);
% b(4).CData = cmap(1,:);
% b(5).CData = cmap(7,:);
% b(6).CData = cmap(6,:);
% b(7).CData = cmap(5,:);
% 
% 
% 

tickFont = 24;
% legendFont = 20;
legendFont = 20;
axisFont = 28;

l=legend('FIXa in solution','FIXa on surface',...
    'FIXa:FX on surface',...
     'FIXa:bsAb:FXa on surface',...
     'FIXa:bsAb:FX on surface',...%'IXa:M:FX on surface',...
     'FIXa:bsAb on surface',...
     'FIXa:bsAb in solution','Location','eastoutside','FontSize',legendFont);

 set(gca,'FontSize',20)
ylim([0 1.1])
ylabel('Fraction of Total FIXa','Fontsize',axisFont)
hold on

h=plot(1,1.05,'.','color',color4,MarkerSize=marker)
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

h=plot(2,1.05,'.','color',color3,MarkerSize=marker)
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

h=plot(3,1.05,'.','color',color2,MarkerSize=marker)
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

h=plot(4,1.05,'.','color',color1,MarkerSize=marker)
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

 set(gca,'FontSize',20)
text(0.025,0.95,chars{2},'Units','normalized','Color',[0,0,0],'FontSize',labelFont,'FontWeight','bold')

 %%%%%%%%%%%%%%%%%%%%%%
ax3 = gca(fig2);
newAx3 = nexttile(t);
text(0.025,0.95,chars{3},'Units','normalized','Color',[1,1,1],'FontSize',labelFont,'FontWeight','bold')

copyAxesProperties(ax3, newAx3);
hold on
kcat=0.485;
plot(kcat,Kd_MX,'.','color',color1,MarkerSize=marker)
colormap(parula)
caxis([0 0.125])
ylabel('K_D^{X,bsAb} (nM)')
set(gca,'FontSize',20)
yticks([1 100 10000]);
yticklabels({'10^0','10^2','10^4'})

ax4 = gca(fig3);
newAx4 = nexttile(t);
text(0.025,0.95,chars{4},'Units','normalized','Color',[1,1,1],'FontSize',labelFont,'FontWeight','bold')

copyAxesProperties(ax4, newAx4);
hold on
%plot(0.75,Kd_IXaM,'.','color',color1,MarkerSize=marker)
colormap(parula)
kcat=0.485;
plot(kcat,Kd_IXaM,'.','color',color1,MarkerSize=marker)
plot(kcat,Kd_IXaM/10,'.','color',color2,MarkerSize=marker)
plot(kcat,Kd_IXaM/100,'.','color',color3,MarkerSize=marker)
plot(kcat,Kd_IXaM/1000,'.','color',color4,MarkerSize=marker)
cb = colorbar
caxis([0 0.125])
cb.Label.String = 'Velocity (nM/s)';
ylabel('K_D^{IXa,bsAb} (nM)')
yticks([1 100 10000]);
yticklabels({'10^0','10^2','10^4'})

set(gca,'FontSize',20)
close(fig1);
close(fig2);
close(fig3);


function copyAxesProperties(sourceAx, targetAx)
    % Copy children
    copyobj(allchild(sourceAx), targetAx);
    % Copy axis properties
    properties = {'XLabel', 'YLabel', 'Title', 'XLim', 'YLim', 'XTick', 'YTick', 'XTickLabel', 'YTickLabel', 'XScale', 'YScale'};
    for i = 1:length(properties)
        targetAx.(properties{i}) = sourceAx.(properties{i});
    end
    % Copy additional properties
    targetAx.XLabel.FontSize = sourceAx.XLabel.FontSize;
    targetAx.YLabel.FontSize = sourceAx.YLabel.FontSize;
    targetAx.Title.FontSize = sourceAx.Title.FontSize;
    targetAx.XLabel.FontWeight = sourceAx.XLabel.FontWeight;
    targetAx.YLabel.FontWeight = sourceAx.YLabel.FontWeight;
    targetAx.Title.FontWeight = sourceAx.Title.FontWeight;
    targetAx.XAxis.FontSize = sourceAx.XAxis.FontSize;
    targetAx.YAxis.FontSize = sourceAx.YAxis.FontSize;
    targetAx.XAxis.FontWeight = sourceAx.XAxis.FontWeight;
    targetAx.YAxis.FontWeight = sourceAx.YAxis.FontWeight;
end