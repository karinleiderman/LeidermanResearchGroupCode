%% emi_plotGen.m
% genrates all matlab figures for the emi paper
% reads in mat files and generates plots in the order they appear in the
% paper. 
% Created: 9/30/24; updated 12/10/24

clear all
close all

%% import data files
% velocity: w/+w/o lipid, FX act by IXa, FX act by TF:VIIa for varying M
vels = importdata("emiVelocities.mat");
vel_7noM = importdata("velFit_noM.mat");
vel_7M = importdata("velFit_M.mat");

% Lipid dependence: IXa, X, M, w/ + w/o lipid
data_fig4 = readmatrix('20240812_emicizumab_c_s_lipid.xlsx','Sheet', 'Sheet1','Range','A24:C35');

% TF:VIIa study
% Varied FX: TF:VIIa, varied FX, lipid, w/ + w/o emi
data_fig5{1} = importdata('oneArm_varyFX.mat');
% Varied M: TF:VIIa, FX, lipid, varied M
data_fig5{2} =  importdata('oneArmFwd.mat');
resMCMC = importdata('oneArmMCMC.mat');
% FIXME!!! break this up or re-run and save to have the same structure as
% the other modelSol
data_fig5{2}.modSol{4}.Ptot = resMCMC;

% IXa study
% no M
% data_fig6{1} = importdata('M=0_2p50000sim_MCMC_FIXaMod_07-30-2024_17-29.mat');
data_fig6{1} = importdata('FIXaMod_noM_fwd.mat');
% data_fig6{2} = importdata('3p100000sim_MCMC_FIXaMod_10-01-2024_06-08.mat');
% data_fig6{2} = importdata('3p50000sim_MCMC_FIXaMod_10-16-2024_14-56.mat');
data_fig6{2} = importdata('FIXaMod_varyM_fwd2.mat');
% validation data: vary FX
data_fig7{1} = importdata('FIXaMod_varyFX_fwd2.mat');

data_fig8{1} = importdata('twoArmfwd_varyKds.mat');
data_fig8{2} = importdata('twoArmfwd_varykcatM_KdMX.mat');
data_fig8{3} = importdata('twoArmfwd_varykcatM_KdIXaM.mat');

% Supplemental figures
data_sup1 = readmatrix('20250321_data.xlsx','Sheet', 'Sheet1','Range','A2:C25');


%% Plot features
tickFont = 34;
legendFont = 20;
axisFont = 48;
marker = 55;
lineWid = 5;

%Colors for plotting
colors = colormap(lines(24));
% cmap = colormap(jet(256));
% % cmap = cmap(200:256,:)
% % non_neon_indices = 50:200;
% 
% % Avoid bright colors (low intensity), focus on deep orange and red
% is_deep_tone = (sum(cmap, 2) < 1.5);
% % Select the indices corresponding to deep tones
% filtered_colors = cmap(is_deep_tone, :);
% 
% 
% % Define criteria for neon green
% is_neon_green = (cmap(:, 1) < 0.2) & (cmap(:, 2) > 0.8) & (cmap(:, 3) < 0.2);
% 
% % Define criteria for neon yellow
% is_neon_yellow = (cmap(:, 1) > 0.8) & (cmap(:, 2) > 0.8) & (cmap(:, 3) < 0.2);
% 
% % valid_indices = find(~is_neon_yellow & ~is_neon_green);
% % valid_indices = 150:250;
% num_colors = 24;
% step = floor(length(valid_indices) / num_colors);
% selected_indices = valid_indices(1:step:end);
% 
% % colors = cmap(selected_indices, :);
% colors = filtered_colors(1:5:end,:);

alphabet = ('A':'Z').';
chars = num2cell(alphabet)';

%% Figure 4: Compare velocities 
tickFont = 28;
legendFont = 20;
axisFont = 24;
marker = 55;
lineWid = 5;
plotCount = 1;
labelFont = 24;

vel_L = vels.vel_L;
MvalsL = vels.MvalsL;
vel_L(8)=[];
MvalsL(8)=[];

common_yticks = linspace(0, 0.18, 5); % Define the same number of ticks
common_xticks = logspace(0, 4, 2); % Define the same number of ticks
figure('Units','inches')
pos = get(gcf,'pos');
set(gcf,'pos',[pos(1) pos(2) 20 6]) % [x y width height]
tile1=tiledlayout(1,3)
nexttile
plot(MvalsL,vel_L,'.','color',colors(2,:),MarkerSize=marker,DisplayName='plus lipid, 0.2nM IXa')
hold on
% plot(vels.MvalsL,vels.vel_L./VmaxL,'-','color',colors(1,:),LineWidth=lineWid,HandleVisibility='off')
plot(vels.MvalsL,vels.vel_noL,'.','color',colors(3,:),MarkerSize=marker,DisplayName='no lipid, 10nM IXa')
% plot(MvalsL,vel_noL./VmaxL,'-','color',colors(2,:),LineWidth=lineWid,HandleVisibility='off')
text(0.025,0.95,chars{plotCount},'Units','normalized','FontSize',labelFont,'FontWeight','bold')
xlim([0 inf])
set(gca,'FontSize',tickFont,'XScale','log','XTick',common_xticks)%,'YTick',common_yticks)
legend('Location','best','FontSize',legendFont)
plotCount=plotCount+1;

nexttile
plot(vels.Mvals9a,vels.vel9a(1:2:end),'.','color',colors(1,:),MarkerSize=marker,DisplayName='1nM IXa')
hold on
% plot(vels.Mvals9a,vels.vel9a(1:2:end),'-','color',colors(4,:),LineWidth=lineWid,HandleVisibility='off')
text(0.025,0.95,chars{plotCount},'Units','normalized','FontSize',labelFont,'FontWeight','bold')
xlim([0 inf])
set(gca,'FontSize',tickFont,'XScale','log','XTick',common_xticks,'YTick',[0,0.05,0.1])
legend('Location','southeast','FontSize',legendFont)
plotCount=plotCount+1;

nexttile
plot(vels.Mvals7a,vels.vel7a,'.','color',colors(5,:),MarkerSize=marker,DisplayName='0.5nM TF:VIIa')
hold on
% plot(vels.Mvals7a,vels.vel7a,'-','color',colors(6,:),LineWidth=lineWid,HandleVisibility='off')
text(0.025,0.95,chars{plotCount},'Units','normalized','FontSize',labelFont,'FontWeight','bold')
xlim([0 inf])
% ylim([0.08 0.18])
ylim([0 0.18])
set(gca,'FontSize',tickFont,'XScale','log','XTick',[100,10000])%,'YTick',[0.08,0.13,0.18])
% sgtitle('No M','FontSize',axisFont)
ylabel(tile1,'Rate of FX activation (nM/s)','FontSize',axisFont)
xlabel(tile1,'Emicizumab (nM)','FontSize',axisFont)
legend('Location','southwest','FontSize',legendFont)

%% Figure 5: TF:VIIa study

% 3 figures of: (1 Full TF:VIIa + M model, (2) model w/o lipid binding, 
% (3) model w/o activation
% M0=50000,3125,390,0 (idx=1,5,8,12)
%Set font size
tickFont = 22;
legendFont = 28;
axisFont = 38;
marker = 55;
lineWid = 5;
labelFont = 28;
plotCount=1;

figure('Units','inches')
pos = get(gcf,'pos');
set(gcf,'pos',[pos(1) pos(2) 13 11]) % [x y width height]

fig5 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
plotCount = 1;

legendHandles = []; % Store only one line per M value for the legend
legendTexts = {};

% Track if we've already added a legend line for each Mval
usedM = containers.Map('KeyType', 'double', 'ValueType', 'logical');

for modIter = 1:4
    colorCount = 4;
    nexttile
    % for i = [12, 8, 5, 1]
    for i = 1:length(data_fig5{2}.Mvals)
        Mval = round(data_fig5{2}.Mvals(i)/1000, 2);
        txt = ['M = ', num2str(Mval), ' \muM'];

        % Dots (experimental data)
        for j = 2:width(data_fig5{2}.data{i}.ydata)
            plot(data_fig5{2}.data{i}.ydata(:,1)./60, ...
                 data_fig5{2}.data{i}.ydata(:,j), '.', ...
                 'Color', colors(colorCount,:), ...
                 'MarkerSize', marker, ...
                 'HandleVisibility', 'off');
            hold on;
        end

        % Line (model data)
        h = plot(data_fig5{2}.data{i}.ydata(:,1)./60, ...
                 data_fig5{2}.modSol{modIter}.Ptot(:,i), '-', ...
                 'Color', colors(colorCount,:), ...
                 'LineWidth', lineWid);

        % Only collect one line per unique Mval for the legend
        if ~isKey(usedM, Mval)
            legendHandles(end+1) = h;
            legendTexts{end+1} = txt;
            usedM(Mval) = true;
            h.DisplayName = txt; % Safe to assign now
        else
            h.HandleVisibility = 'off'; % Hide in legend
        end

        colorCount = colorCount + 1;
    end

    text(0.025, 0.95, chars{plotCount}, ...
         'Units', 'normalized', ...
         'FontSize', labelFont, ...
         'FontWeight', 'bold');
    plotCount = plotCount + 1;

    if mod(modIter, 2) == 0
        set(gca, 'FontSize', tickFont, 'YTick', []);
    else
        set(gca, 'FontSize', tickFont);
    end
end

xlabel(fig5, 'Time (min)', 'FontSize', axisFont);
ylabel(fig5, 'FXa (nM)', 'FontSize', axisFont);

% lgd = legend(legendHandles, legendTexts, ...
%     'Location', 'northoutside', ...
%     'Orientation', 'horizontal', ...
%     'FontSize', legendFont, ...
%     'NumColumns', numel(legendHandles));
lgd = legend(legendHandles, legendTexts, ...
    'Location', 'eastoutside', ...
    'Orientation', 'horizontal');

% Attach the legend to the top of the tiled layout
lgd.Layout.Tile = 'north';


%% Figure 6: IXa Study - parameter estimation results
% 2 figures of: (1) Parameter estimation results for no M, (2) results with
% M
% M0=50000,3125,390,0 (idx=1,5,8,12)
%Set font size
plotCount = 1;
tickFont = 28; % tickmarks, i.e. numbers on axis
legendFont = 20;
axisFont = 28; % axis label, e.g. Time (min)
marker = 25;
lineWid = 2;
labelFont = 28; % plot label, i.e. A, B, C
plotCount=1;
colorCount=1;

figure('Units','inches')
pos = get(gcf,'pos');
set(gcf,'pos',[pos(1) pos(2) 20 7]) % [x y width height]
fig6 = tiledlayout(1,2);
nexttile
% for j=1:length(data_fig6{1}.data)
for j=1:length(data_fig6{1}.data)
    Xval = data_fig6{1}.data{j}.y0(1) + data_fig6{1}.data{j}.y0(9) + data_fig6{1}.data{j}.y0(15);
    txt = ['X=',num2str(floor(Xval)),'nM'];
    plot(data_fig6{1}.data{j}.ydata(1:10:end,1)./60,data_fig6{1}.data{j}.ydata(1:10:end,2),'.','Color',colors(j,:),MarkerSize=marker,HandleVisibility='off')
    hold on
    plot(data_fig6{1}.data{j}.ydata(1:10:end,1)./60,data_fig6{1}.data{j}.ydata(1:10:end,3),'.','Color',colors(j,:),MarkerSize=marker,HandleVisibility='off')
    plot(data_fig6{1}.data{j}.ydata(:,1)./60,data_fig6{1}.C(:,j),'-','Color',colors(j,:),LineWidth=lineWid,DisplayName=txt)
end
% xlim([0 95])
ylim([-inf 4.1*10^5])
set(gca,'FontSize',tickFont)
legend('Location','bestoutside','FontSize',legendFont)
text(0.025,0.95,chars{plotCount},'Units','normalized','FontSize',labelFont,'FontWeight','bold')
plotCount=plotCount+1;
nexttile
% for i=1:2:length(data_fig6{2}.data)
for i=[1,5,13,15,17,19,21]
    Mval = data_fig6{2}.data{i}.y0(12) + data_fig6{2}.data{i}.y0(13)+ data_fig6{2}.data{i}.y0(15);
    txt = ['M=',num2str(round(Mval,2)),'nM'];
    plot(data_fig6{2}.data{i}.ydata(1:4:end,1)./60,data_fig6{2}.data{i}.ydata(1:4:end,2),'.','Color',colors(i,:),MarkerSize=marker,HandleVisibility='off')
    hold on
    plot(data_fig6{2}.data{i+1}.ydata(1:4:end,1)./60,data_fig6{2}.data{i+1}.ydata(1:4:end,2),'.','Color',colors(i,:),MarkerSize=marker,HandleVisibility='off')
    plot(data_fig6{2}.data{i}.ydata(:,1)./60,data_fig6{2}.C(:,i),'-','Color',colors(i,:),LineWidth=lineWid,DisplayName=txt)
end
ylim([-inf 4.1*10^5])
set(gca,'FontSize',tickFont)
legend('Location','bestoutside','FontSize',legendFont)
text(0.025,0.95,chars{plotCount},'Units','normalized','FontSize',labelFont,'FontWeight','bold')
xlabel(fig6,'Time (min)','fontsize',axisFont)
ylabel(fig6,'Chromogenic Substrate (nM)','fontsize',axisFont)

%% Figure 7: IXa Study - validation with addtional data
% 1 tiled figure of model fit to data for varying FX
tickFont = 22; % tickmarks, i.e. numbers on axis
legendFont = 12;
axisFont = 28; % axis label, e.g. Time (min)
marker = 35;
lineWid = 3;
labelFont = 28;
plotCount=1;
colorCount=1;


% common_xticks = linspace(0, 10, 3); % Define the same number of ticks
common_xticks = [0,2,4,6,8]; % Define the same number of ticks

figure('Units','inches')
pos = get(gcf,'pos');
set(gcf,'pos',[pos(1) pos(2) 20 8]) % [x y width height]
tile1 = tiledlayout(2,4,"TileSpacing","tight");
% tile1 = tiledlayout(2,4);
for i=1:2:length(data_fig7{1}.data)-4
    nexttile 
    Xval = data_fig7{1}.data{i}.y0(1) + data_fig7{1}.data{i}.y0(9) + data_fig7{1}.data{i}.y0(15);
    % Mval = data_fig7{1}.data{i}.y0(12) + data_fig7{1}.data{i}.y0(13)+ data_fig7{1}.data{i}.y0(15);
    txt = ['FX = ',num2str(floor(Xval)),'nM'];
    plot(data_fig7{1}.data{i}.ydata(1:48,1)./60,data_fig7{1}.data{i}.ydata(1:48,2)./1000,'.','Color',colors(i,:),'MarkerSize',10,'HandleVisibility','off')%'DisplayName',txt)
    hold on
    plot(data_fig7{1}.data{i+1}.ydata(1:48,1)./60,data_fig7{1}.data{i+1}.ydata(1:48,2)./1000,'.','Color',colors(i,:),'MarkerSize',10,'HandleVisibility','off')%'DisplayName',txt)
    plot(data_fig7{1}.data{i}.ydata(1:48,1)./60,data_fig7{1}.C(1:48,i)./1000,'-','Color',colors(i,:),'LineWidth',lineWid,'DisplayName','Model')
    text(4,100,txt,'fontsize',18,'FontWeight','normal')
    ylim([0 400]);
    xlim([0 11]);
    if (plotCount~=1) & (plotCount~=5) 
        set(gca,'FontSize',tickFont,'XTick',common_xticks,'ytick',[])
    end
    if plotCount<5
        set(gca,'FontSize',tickFont,'xtick',[])
    else
        set(gca,'FontSize',tickFont,'XTick',common_xticks)
    end
    plotCount = plotCount + 1;
end
% set(gca,'FontSize',tickFont)
% xlim([tdata(1)./60-1 inf]);
% legend('Location','bestoutside','FontSize',legendFont)
% title('Data and Model','FontSize',axisFont)
xlabel(tile1,'Time (min)','FontSize',axisFont)
ylabel(tile1,'Chromogenic Substrate (\muM)','FontSize',axisFont)

%% Figure 8: IXa Study - model analysis 
% % tiled layout of 3 heatmaps of velocity with varying parameters
% plotCount=1;
% Kd_IXaM=5.5276e+03;
% Kd_MX=55.7235;
% kcatEST=0.75;
% marker = 55;
% 
% common_ticks = logspace(0, 4, 3); % Define the same number of ticks
% 
% figure('Units','inches')
% pos = get(gcf,'pos');
% set(gcf,'pos',[pos(1) pos(2) 20 7]) % [x y width height]
% tile8 = tiledlayout(1,3,'TileSpacing','tight');
% % KD_IXaM vs. KD_MX
% nexttile
% pcolor(data_fig8{1}.varValues(1,:),data_fig8{1}.varValues(2,:),data_fig8{1}.V')
% hold on
% shading flat;
% caxis([0,0.18]);
% plot(Kd_IXaM,Kd_MX,'r.',MarkerSize=marker)
% plot(data_fig8{1}.varValues(1,:),data_fig8{1}.varValues(2,:),'k-',LineWidth=2)
% 
% text(0.025,0.95,chars{plotCount},'Units','normalized','Color',[1,1,1],'FontSize',labelFont,'FontWeight','bold')
% 
% set(gca,'FontSize',tickFont,'XScale','log','YScale','log','XTick',common_ticks,'YTick',common_ticks)
% xlabel('K_D^{IXa,M} (nM)','fontsize',axisFont)
% ylabel('K_D^{M,X} (nM)','fontsize',axisFont)
% plotCount = plotCount + 1;
% % kcatM vs. KD_MX
% nexttile
% pcolor(data_fig8{2}.varValues(1,:),data_fig8{2}.varValues(2,:),data_fig8{2}.V')
% hold on
% shading flat;
% caxis([0,0.18]);
% plot(kcatEST,Kd_MX,'r.',MarkerSize=marker)
% 
% text(0.025,0.95,chars{plotCount},'Units','normalized','Color',[1,1,1],'FontSize',labelFont,'FontWeight','bold')
% 
% set(gca,'FontSize',tickFont,'XScale','log','YScale','log')
% ylabel('K_D^{M,X} (nM)','fontsize',axisFont)
% xlabel('k^{cat}(s^{-1})','fontsize',axisFont)
% plotCount = plotCount + 1;
% % kcatM vs. KD_IXaM
% nexttile
% pcolor(data_fig8{3}.varValues(1,:),data_fig8{3}.varValues(2,:),data_fig8{3}.V')
% hold on
% shading flat;
% plot(kcatEST,Kd_IXaM,'r.',MarkerSize=marker)
% 
% text(0.025,0.95,chars{plotCount},'Units','normalized','Color',[1,1,1],'FontSize',labelFont,'FontWeight','bold')
% 
% set(gca,'FontSize',tickFont,'XScale','log','YScale','log')
% ylabel('K_D^{IXa,M} (nM)','fontsize',axisFont)
% xlabel('k^{cat}(s^{-1})','fontsize',axisFont)
% cb = colorbar;
% cb.Layout.Tile = 'east';
% cb.Ticks = [0, 0.06, 0.12, 0.18];
% colormap parula
% caxis([0,0.18]);
%% Figure 8: IXa Study - model analysis (individual plots)
compilefigs


% tiled layout of 3 heatmaps of velocity with varying parameters
% plotCount=1;
% Kd_IXaM=5.5276e+03;
% Kd_MX=55.7235;
% kcatEST=0.75;
% marker = 55;
% 
% common_ticks = logspace(0, 4, 3); % Define the same number of ticks
% 
% figure('Units','inches')
% pos = get(gcf,'pos');
% set(gcf,'pos',[pos(1) pos(2) 10 7]) % [x y width height]
% % tile8 = tiledlayout(1,3,'TileSpacing','tight');
% % KD_IXaM vs. KD_MX
% % nexttile
% pcolor(data_fig8{1}.varValues(1,:),data_fig8{1}.varValues(2,:),data_fig8{1}.V')
% hold on
% shading flat;
% cb = colorbar;
% % cb.Layout.Tile = 'east';
% cb.Ticks = [0, 0.06, 0.12, 0.18];
% colormap parula
% caxis([0,0.18]);
% plot(Kd_IXaM,Kd_MX,'r.',MarkerSize=marker)
% plot(data_fig8{1}.varValues(1,:),data_fig8{1}.varValues(2,:),'k-',LineWidth=2)
% 
% text(0.025,0.95,chars{plotCount},'Units','normalized','Color',[1,1,1],'FontSize',labelFont,'FontWeight','bold')
% 
% set(gca,'FontSize',tickFont,'XScale','log','YScale','log','XTick',common_ticks,'YTick',common_ticks)
% xlabel('K_D^{IXa,M} (nM)','fontsize',axisFont)
% ylabel('K_D^{M,X} (nM)','fontsize',axisFont)
% plotCount = plotCount + 1;
% % kcatM vs. KD_MX
% figure('Units','inches')
% pos = get(gcf,'pos');
% set(gcf,'pos',[pos(1) pos(2) 10 7]) % [x y width height]
% pcolor(data_fig8{2}.varValues(1,:),data_fig8{2}.varValues(2,:),data_fig8{2}.V')
% hold on
% shading flat;
% cb = colorbar;
% % cb.Layout.Tile = 'east';
% cb.Ticks = [0, 0.06, 0.12, 0.18];
% colormap parula
% caxis([0,0.18]);
% plot(kcatEST,Kd_MX,'r.',MarkerSize=marker)
% 
% % text(0.025,0.95,chars{plotCount},'Units','normalized','Color',[1,1,1],'FontSize',labelFont,'FontWeight','bold')
% 
% set(gca,'FontSize',tickFont,'XScale','log','YScale','log')
% ylabel('K_D^{M,X} (nM)','fontsize',axisFont)
% xlabel('k^{cat}(s^{-1})','fontsize',axisFont)
% plotCount = plotCount + 1;
% % kcatM vs. KD_IXaM
% figure('Units','inches')
% pos = get(gcf,'pos');
% set(gcf,'pos',[pos(1) pos(2) 10 7]) % [x y width height]
% pcolor(data_fig8{3}.varValues(1,:),data_fig8{3}.varValues(2,:),data_fig8{3}.V')
% hold on
% shading flat;
% plot(kcatEST,Kd_IXaM,'r.',MarkerSize=marker)
% 
% % text(0.025,0.95,chars{plotCount},'Units','normalized','Color',[1,1,1],'FontSize',labelFont,'FontWeight','bold')
% 
% set(gca,'FontSize',tickFont,'XScale','log','YScale','log')
% ylabel('K_D^{IXa,M} (nM)','fontsize',axisFont)
% xlabel('k^{cat}(s^{-1})','fontsize',axisFont)
% cb = colorbar;
% % cb.Layout.Tile = 'east';
% cb.Ticks = [0, 0.06, 0.12, 0.18];
% colormap parula
% caxis([0,0.18]);

%% Supplemental Figures
tickFont = 28;
legendFont = 28;
axisFont = 24;
marker = 55;
lineWid = 5;
plotCount = 1;
labelFont = 24;

figure('Units','inches')
pos = get(gcf,'pos');
set(gcf,'pos',[pos(1) pos(2) 13 11]) % [x y width height]
plot(data_sup1(:,1),data_sup1(:,2),'.','color',colors(2,:),MarkerSize=marker,DisplayName='FXa and M')
hold on
plot(data_sup1(:,1),data_sup1(:,3),'.','color',colors(1,:),MarkerSize=marker,DisplayName='FXa')
set(gca,'FontSize',tickFont,'XScale','log')%,'YScale','log')
ylabel('Fraction FXa in Solution','fontsize',axisFont)
xlabel('Lipid (\muM)','fontsize',axisFont)
legend('Location','best','FontSize',legendFont)

%% EXTRA PLOTS
% velocity for TF:VIIa+FX w + w/o M
% figure('Units','inches')
% pos = get(gcf,'pos');
% set(gcf,'pos',[pos(1) pos(2) 20 11]) % [x y width height]
% plot(data_fig4(:,1),data_fig4(:,2),'k.','MarkerSize',marker,'DisplayName','M = 0 \mu M');
% hold on
% plot(data_fig4(:,1),data_fig4(:,3),'b.','MarkerSize',marker,'DisplayName','M = 50 \mu M');
% set(gca,'FontSize',tickFont)
% xlabel('Time (s)','fontsize',axisFont)
% ylabel('Rate of FXa generation (nM/min)','fontsize',axisFont)
% text(0.025,0.95,chars{plotCount},'Units','normalized','FontSize',labelFont,'FontWeight','bold')

% FIXa Mod plots with MCMC results (not fwd run)
% figure('Units','inches')
% pos = get(gcf,'pos');
% set(gcf,'pos',[pos(1) pos(2) 20 11]) % [x y width height]
% fig6 = tiledlayout(1,2);
% nexttile
% for j=3:length(data_fig6{1}.data)
%     Xval = data_fig6{1}.data{j}.y0(1) + data_fig6{1}.data{j}.y0(9) + data_fig6{1}.data{j}.y0(15);
%     txt = ['X=',num2str(floor(Xval)),'nM'];
%     plot(data_fig6{1}.data{j}.ydata(:,1)./60,data_fig6{1}.data{j}.ydata(:,2),'.','Color',colors(j,:),MarkerSize=marker,HandleVisibility='off')
%     hold on
%     plot(data_fig6{1}.data{j}.ydata(:,1)./60,data_fig6{1}.C(:,j),'-','Color',colors(j,:),LineWidth=lineWid,DisplayName=txt)
% end
% set(gca,'FontSize',tickFont)
% legend('Location','bestoutside','FontSize',legendFont)
% text(0.025,0.95,chars{plotCount},'Units','normalized','FontSize',labelFont,'FontWeight','bold')
% plotCount=plotCount+1;
% nexttile
% for i=1:2:length(data_fig6{2}.data)
%     Mval = data_fig6{2}.data{i}.y0(12) + data_fig6{2}.data{i}.y0(13)+ data_fig6{2}.data{i}.y0(15);
%     txt = ['M=',num2str(round(Mval,2)),'nM'];
%     plot(data_fig6{2}.data{i}.ydata(:,1)./60,data_fig6{2}.data{i}.ydata(:,2),'.','Color',colors(i,:),MarkerSize=marker,HandleVisibility='off')
%     hold on
%     plot(data_fig6{2}.data{i+1}.ydata(:,1)./60,data_fig6{2}.data{i+1}.ydata(:,2),'.','Color',colors(i,:),MarkerSize=marker,HandleVisibility='off')
%     plot(data_fig6{2}.data{i}.ydata(:,1)./60,data_fig6{2}.sol{i}.MCMC,'-','Color',colors(i,:),LineWidth=lineWid,DisplayName=txt)
% end
% set(gca,'FontSize',tickFont)
% legend('Location','bestoutside','FontSize',legendFont)
% text(0.025,0.95,chars{plotCount},'Units','normalized','FontSize',labelFont,'FontWeight','bold')
% xlabel(fig6,'Time (min)','fontsize',axisFont)
% ylabel(fig6,'C (nM)','fontsize',axisFont)

% for i=1:length(data)
%     % nexttile 
%     Xval = data{i}.y0(1) + data{i}.y0(9) + data{i}.y0(15);
%     Mval = data{i}.y0(12) + data{i}.y0(13)+ data{i}.y0(15);
%     txt = ['M=',num2str(round(Mval,2)),'nM;','X=',num2str(floor(Xval)),'nM'];
%     plot(data{i}.ydata(:,1)./60,data{i}.ydata(:,2),'.','Color',colors(i,:),'MarkerSize',10,'HandleVisibility','off')%'DisplayName',txt)
%     hold on
%     % plot(data{i+1}.ydata(:,1)./60,data{i+1}.ydata(:,2),'.','Color',colors(i,:),'MarkerSize',10,'HandleVisibility','off')%'DisplayName',txt)
%     plot(data{i}.ydata(:,1)./60,C(:,i),'-','Color',colors(i,:),'LineWidth',lineWid,'DisplayName','Model')
%     % plot(data{i+1}.ydata(:,1)./60,C(:,i+1),'-','Color',colors(i,:),'LineWidth',lineWid,'DisplayName','Model')
%     % ylim([0 4.1*10^5]);
%     % xlim([0 10]);
% end
% % set(gca,'FontSize',tickFont)
% % xlim([tdata(1)./60-1 inf]);
% % legend('Location','bestoutside','FontSize',legendFont)
% % title('Data and Model','FontSize',axisFont)
% xlabel(tile1,'Time (min)','FontSize',axisFont)
% ylabel(tile1,'C (nM)','FontSize',axisFont)