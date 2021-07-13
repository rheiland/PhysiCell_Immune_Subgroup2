load('simulated_data.mat')
timetotal = 72;
hours_to_days = linspace(0,15,timetotal);
 
 figure
 hold on 
 plot(hours_to_days,uninfected,'LineWidth',2)
 plot(hours_to_days,infected,'--','LineWidth',2)
 plot(hours_to_days,dead,':','LineWidth',2)
 plot(hours_to_days,antiviral_cells,':','LineWidth',2)
 ylabel('Number of cells')
% ylim([0 3000])
 yyaxis right
 plot(hours_to_days,virion','-.','LineWidth',2)
 ylabel('Total virions')
 legend('Uninfected cells','Infected cells','Dead cells','Antiviral cells','Virions')
 xlabel('Time (days)')
 set(gca,'FontSize',15)
 ax = gca;
 ax.YAxis(1).Color = 'k';
 ax.YAxis(2).Color = 'k';
 saveas(gcf,'F1.png')

 figure
 hold on
 plot(hours_to_days,macrophagesinactive,'Color',[0.04 0.31 0.49],'LineWidth',2)
 plot(hours_to_days,macrophagesactive,'--','Color',[.3 .75 .93],'LineWidth',2)
 plot(hours_to_days,neutrophils,':','Color',[.47 .67 .19],'LineWidth',2)
 legend('Macrophages (inactive)','Macrophages (active)','Neutrophils')
 xlabel('Time (days)')
 ylabel('Number of cells')
 set(gca,'FontSize',15)
 saveas(gcf,'F2.png')
 
 
 figure
 hold on
 plot(hours_to_days,macrophagesinactive+macrophagesactive,'Color',[0.04 0.31 0.49],'LineWidth',2)
 plot(hours_to_days,neutrophils,':','Color',[.47 .67 .19],'LineWidth',2)
 plot(hours_to_days,DC,'Color',[1 0 0],'LineWidth',2)
 legend('Total Macrophages','Neutrophils','Dendritic cells')
 xlabel('Time (days)')
 ylabel('Number of cells')
 set(gca,'FontSize',15)
 saveas(gcf,'F2.png')
 
figure
 hold on
 plot(hours_to_days,CD4,'--','Color',[1 0.07 0.65],'LineWidth',2)
 plot(hours_to_days,CD8,':','Color',[.64 .08 .18],'LineWidth',2)
  legend('DCs','CD8 T cells','CD4 T cells')
 xlabel('Time (days)')
 ylabel('Number of cells')
 set(gca,'FontSize',15)
 saveas(gcf,'F3.png')
 
 figure
 hold on 
 plot(hours_to_days,cytokine,'Color',[1 0.41 0.16],'LineWidth',2)
 plot(hours_to_days,chemokine,'--','Color',[0.12 0.64 0.54],'LineWidth',2)
 plot(hours_to_days,debris,':','Color',[.73 0.46 0.9],'LineWidth',2)
 plot(hours_to_days,IFN,':','LineWidth',2)
 plot(hours_to_days,ROS,':','LineWidth',2)
 ylabel('Total substrates')
 legend('Pro-inflammatory cytokine','Chemokine','Debris','IFN','ROS')
 xlabel('Time (days)')
 set(gca,'FontSize',15)
 saveas(gcf,'F4.png')
 
 
 %% load Amber's data

 load('AmberData.mat')
 
 domain_ml = 0.0008;
 
 figure
 hold on 
 yyaxis left
 plot(time_virus,10.^ASvirus_data,'o:','Color',[135,135,135]/255,'LineWidth',1)
ylabel('Virion (log_{10}(copies/ml))')
 yyaxis right
 plot(hours_to_days,virion'/domain_ml,'-.','Color',[214,96,77]/255,'LineWidth',2)
 xlim([0.1 12])
ylabel('Virion (copies/ml)')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('Virion')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [214,96,77]/255;

 figure
 hold on 
 yyaxis left
 plot(time_virus,10.^ASvirus_data,'o:','Color',[135,135,135]/255,'LineWidth',1)
set(gca,'yscale','log')
ylabel('Virion (copies/ml)')
 yyaxis right
 plot(hours_to_days,virion'/domain_ml,'-.','Color',[214,96,77]/255,'LineWidth',2)
set(gca,'yscale','log')
 xlim([0.1 12])
ylabel('Virion (copies/ml)')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('Virion')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [214,96,77]/255;
 
figure
hold on 
yyaxis left
plot(time_cells,ASneutrophils,'o:','Color',[135,135,135]/255,'LineWidth',1)
 ylabel('Cells')
yyaxis right
plot(hours_to_days,neutrophils,'-.','Color',[67,147,195]/255,'LineWidth',2)
 ylabel('Cells')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('Neutrophils')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [67,147,195]/255;

figure
hold on 
yyaxis left
plot(time_cells,ASmacs_tot,'o:','Color',[135,135,135]/255,'LineWidth',1)
 ylabel('Cells')
yyaxis right
plot(hours_to_days,macrophagesinactive+macrophagesactive,'-.','Color',[253,174,97]/255,'LineWidth',2)
 ylabel('Cells')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('Macrophages')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [253,174,97]/255;
 
figure
hold on 
yyaxis left
plot(time_cells,ASDCs_tot,'o:','Color',[135,135,135]/255,'LineWidth',1)
 ylabel('Cells')
yyaxis right
plot(hours_to_days,DC,'-.','Color',[90,174,97]/255,'LineWidth',2)
 ylabel('Cells')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('DCs')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [90,174,97]/255;
 
figure
hold on 
yyaxis left
plot(time_cells,ASCD4s,'o:','Color',[135,135,135]/255,'LineWidth',1)
ylabel('Cells')
yyaxis right
 plot(hours_to_days,CD4,'-.','Color',[1 0.07 0.65],'LineWidth',2)
ylabel('Cells')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('CD4+ T cells')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [1 0.07 0.65];

figure
hold on 
yyaxis left
plot(time_cells,ASCD8s,'o:','Color',[135,135,135]/255,'LineWidth',1)
ylabel('Cells')
yyaxis right
plot(hours_to_days,CD8,'-.','Color',[.64 .08 .18],'LineWidth',2)
ylabel('Cells')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('CD8+ T cells')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [.64 .08 .18];
 

figure
hold on 
yyaxis left
plot(time_cytokine,ASIL6,'o:','Color',[135,135,135]/255,'LineWidth',1)
ylabel('log_{10}')
yyaxis right
plot(hours_to_days,cytokine,'-.','Color',[53,151,143]/255,'LineWidth',2)
ylabel('pg')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('IL-6')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [53,151,143]/255;
 

%% Topanta Ross data comparison

load('Topanta_Ross.mat')


 figure
 hold on 
 yyaxis left
 plot(time_virusT,virus_adult,'o:','Color',[135,135,135]/255,'LineWidth',1)
 plot(time_virusT,virus_aged,'o:','Color',[135,135,135]/255,'LineWidth',1)
ylabel('Virion (copies/ml)')
 yyaxis right
 plot(hours_to_days,virion'/domain_ml,'-.','Color',[214,96,77]/255,'LineWidth',2)
 xlim([0.1 12])
ylabel('Virion (copies/ml)')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('Virion')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [214,96,77]/255;

figure
hold on 
yyaxis left
plot(time_macs,macs_adult,'o:','Color',[135,135,135]/255,'LineWidth',1)
plot(time_macs,macs_aged,'o:','Color',[135,135,135]/255,'LineWidth',1)
yyaxis right
plot(hours_to_days,macrophagesactive,'-.','Color',[253,174,97]/255,'LineWidth',2)
ylabel('Cells')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('Macrophages')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [253,174,97]/255;
 
figure
hold on 
yyaxis left
plot(time_macs,DCs_adult,'o:','Color',[135,135,135]/255,'LineWidth',1)
plot(time_macs,DCs_aged,'o:','Color',[135,135,135]/255,'LineWidth',1)
yyaxis right
plot(hours_to_days,DC,'-.','Color',[90,174,97]/255,'LineWidth',2)
 ylabel('Cells')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('DCs')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [90,174,97]/255;
 
figure
hold on 
yyaxis left
plot(time_macs,CD4_adult,'o:','Color',[135,135,135]/255,'LineWidth',1)
plot(time_macs,CD4_aged,'o:','Color',[135,135,135]/255,'LineWidth',1)
yyaxis right
 plot(hours_to_days,CD4,'--','Color',[1 0.07 0.65],'LineWidth',2)
ylabel('Cells')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('CD4+ T cells')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [1 0.07 0.65];


figure
hold on 
yyaxis left
plot(time_macs,CD8_adult,'o:','Color',[135,135,135]/255,'LineWidth',1)
plot(time_macs,CD8_aged,'o:','Color',[135,135,135]/255,'LineWidth',1)
ylabel('Cells')
yyaxis right
plot(hours_to_days,CD8,':','Color',[.64 .08 .18],'LineWidth',2)
ylabel('Cells')
xlabel('Time (days)')
set(gca,'FontSize',15)
title('CD8+ T cells')
ax = gca
ax.YAxis(1).Color = [77,77,77]/255;
ax.YAxis(2).Color = [.64 .08 .18];
 