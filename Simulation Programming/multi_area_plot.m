
clc;clear all;close all
%area 1-markov
load('markov_data003.mat')
%area 2-markov
load('markov_data004.mat')
%multi-area
load('multi-Markov_results.mat')
%area 1
load('area1_Diao_initial_K.mat')
load('area1_mode4_opt_data.mat')
load('area1_pkl_data.mat')
load ('area1_data.mat')
%area 2
load('area2_Diao_initial_K.mat')
load('area2_mode4_opt_data.mat')
load('area2_pkl_data.mat')
load ('area2_data.mat')


figure (1)
subplot(311)
stairs(Jump_time1,sttout1,'LineWidth',1.5)
ax = gca; 
ax.TickDir = 'out';

set(gca,'ytick',[0:1:3]) 
set(gca,'xtick') 
xlabel('(a)  Time Step')
ylabel('Jumping Mode')

title('Area 1');
ylim([0, 3]);
yticks(1:3);

subplot(312)
stairs(Jump_time2,sttout2,'LineWidth',1.5)
ax = gca; 
ax.TickDir = 'out';

set(gca,'ytick',[0:1:3]) 
set(gca,'xtick') 
xlabel('(b)  Time Step')
ylabel('Jumping Mode')

title('Area 2');
ylim([0, 3]);
yticks(1:3);

subplot(313)

stairs(Jump_time,sttout,'LineWidth',1.5)
ax = gca; 
ax.TickDir = 'out';
set(gca,'ytick',[0:1:5]) 
set(gca,'xtick') 
xlabel('(c)  Time Step')
ylabel('Jumping Mode')

title('Joint jumping mode');
ylim([0, 5]);
yticks(1:5);


figure(41)
subplot(221)
g=plot(0:length(a1_ww1_mat)-1,a1_ww1_mat',':')
g(1).Color = 'g'; 
g(5).Color = [0.9290 0.6940 0.1250]; 
g(3).Color = [0 0 1]; 
g(4).Color = [1 0 1]; 
g(2).Color = [0 0.4470 0.7410]; 
 g(1).LineStyle='-'
 g(1).LineWidth = 0.8;
 g(5).LineStyle='--'
 g(5).LineWidth = 0.8;
 g(3).LineStyle=':'
 g(3).LineWidth = 0.8;
 g(4).LineStyle='-.'
 g(4).LineWidth = 0.8;
 g(2).LineStyle='none'
 g(2).Marker='p'
 g(2).LineWidth = 0.8;
 g(2).MarkerSize=1.5
xlabel('Time Steps')
ax = gca; 
ax.TickDir = 'out';
ylim([-120 5])
xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'r--','LineWidth',0.5,'HandleVisibility','off')
ax = gca; 
ax.TickDir = 'out';
title('Mode=1')
hold on

subplot(222)
gg=plot(0:length(a1_ww2_mat)-1,a1_ww2_mat',':')
gg(1).Color = 'g'; 
gg(5).Color = [0.9290 0.6940 0.1250]; 
gg(3).Color = [0 0 1]; 
gg(4).Color = [1 0 1]; 
gg(2).Color = [0 0.4470 0.7410]; 
gg(1).LineStyle='-'
 gg(1).LineWidth = 0.8;
 gg(5).LineStyle='--'
 gg(5).LineWidth = 0.8;
 gg(3).LineStyle=':'
 gg(3).LineWidth = 0.8;
 gg(4).LineStyle='-.'
 gg(4).LineWidth = 0.8;
  gg(2).LineStyle='none'
 gg(2).Marker='p'
 gg(2).LineWidth = 0.8;
 gg(2).MarkerSize=1.5
xlabel(' Time Steps')
ax = gca; 
ax.TickDir = 'out';
ylim([-80 5])
xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'r--','LineWidth',0.5,'HandleVisibility','off')
ax = gca; 

ax.TickDir = 'out';
title('Mode=2')
hold on


subplot(223)
gg=plot(0:length(a1_ww3_mat)-1,a1_ww3_mat',':')
gg(1).Color = 'g';
gg(5).Color = [0.9290 0.6940 0.1250]
gg(3).Color = [0 0 1];
gg(4).Color = [1 0 1]; 
gg(2).Color = [0 0.4470 0.7410]; 
 gg(1).LineStyle='-'
 gg(1).LineWidth = 0.8;
 gg(5).LineStyle='--'
 gg(5).LineWidth = 0.8;
 gg(3).LineStyle=':'
 gg(3).LineWidth = 0.8;
 gg(4).LineStyle='-.'
 gg(4).LineWidth = 0.8;
 gg(2).LineStyle='none'
 gg(2).Marker='p'
 gg(2).LineWidth = 0.8;
 gg(2).MarkerSize=1.5
xlabel(' Time Steps')
ax = gca; 
ax.TickDir = 'out';
ylim([-120 5])
xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'r--','LineWidth',0.5,'HandleVisibility','off')
ax = gca; 
ax.TickDir = 'out';
title('Mode=3')
hold on


subplot(224)
gg=plot(0:length(a1_ww4_mat)-1,a1_ww4_mat',':')
 gg(1).Color = 'g'; 
 gg(5).Color = [0.9290 0.6940 0.1250];
 gg(3).Color = [0 0 1];
 gg(4).Color = [1 0 1];
 gg(2).Color = [0 0.4470 0.7410]; 
 gg(1).LineStyle='-'
 gg(5).LineStyle='--'
 gg(5).LineWidth = 0.8;
 gg(3).LineStyle=':'
 gg(3).LineWidth = 0.8;
 gg(4).LineStyle='-.'
gg(4).LineWidth = 0.8;
  gg(2).LineStyle='none'
 gg(2).Marker='p'
 gg(2).LineWidth = 0.8;
 gg(2).MarkerSize=1.5
xlabel(' Time Steps')
ax = gca; 
ax.TickDir = 'out';
ylim([-85 5])
xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'r--','LineWidth',0.5,'HandleVisibility','off')
ax = gca; 

ax.TickDir = 'out';
legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5','fontsize',6,'location','Best','NumColumns',2)
title('Mode=4')% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
figure(42)

subplot(221)
g=plot(0:length(a2_ww1_mat)-1,a2_ww1_mat',':')
g(1).Color = 'g'; 
g(5).Color = [0.9290 0.6940 0.1250]; 
g(3).Color = [0 0 1]; 
g(4).Color = [1 0 1]; 
g(2).Color = [0 0.4470 0.7410]; 
 g(1).LineStyle='-'
 g(1).LineWidth = 0.8;
 g(5).LineStyle='--'
 g(5).LineWidth = 0.8;
 g(3).LineStyle=':'
 g(3).LineWidth = 0.8;
 g(4).LineStyle='-.'
 g(4).LineWidth = 0.8;
 g(2).LineStyle='none'
 g(2).Marker='p'
 g(2).LineWidth = 0.8;
 g(2).MarkerSize=1.5
xlabel('Time Steps')
ax = gca; 
ax.TickDir = 'out';
ylim([-80 5])
xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'r--','LineWidth',0.5,'HandleVisibility','off')
ax = gca; 

ax.TickDir = 'out';
title('Mode=1')
legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5','fontsize',6,'location','Best','NumColumns',2)


hold on


subplot(222)
gg=plot(0:length(a2_ww2_mat)-1,a2_ww2_mat',':')
gg(1).Color = 'g';
gg(5).Color = [0.9290 0.6940 0.1250]; 
gg(3).Color = [0 0 1]; 
gg(4).Color = [1 0 1];
gg(2).Color = [0 0.4470 0.7410]; 
gg(1).LineStyle='-'
 gg(1).LineWidth = 0.8;
 gg(5).LineStyle='--'
 gg(5).LineWidth = 0.8;
 gg(3).LineStyle=':'
 gg(3).LineWidth = 0.8;
 gg(4).LineStyle='-.'
 gg(4).LineWidth = 0.8;
  gg(2).LineStyle='none'
 gg(2).Marker='p'
 gg(2).LineWidth = 0.8;
 gg(2).MarkerSize=1.5
xlabel(' Time Steps')
ax = gca; 
ax.TickDir = 'out';
ylim([-45 5])
xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'r--','LineWidth',0.5,'HandleVisibility','off')
ax = gca; 

ax.TickDir = 'out';
title('Mode=2')
hold on


subplot(223)
gg=plot(0:length(a2_ww3_mat)-1,a2_ww3_mat',':')
gg(1).Color = 'g'; 
gg(5).Color = [0.9290 0.6940 0.1250]
gg(3).Color = [0 0 1]; 
gg(4).Color = [1 0 1]; 
gg(2).Color = [0 0.4470 0.7410]; 
 gg(1).LineStyle='-'
 gg(1).LineWidth = 0.8;
 gg(5).LineStyle='--'
 gg(5).LineWidth = 0.8;
 gg(3).LineStyle=':'
 gg(3).LineWidth = 0.8;
 gg(4).LineStyle='-.'
 gg(4).LineWidth = 0.8;
 gg(2).LineStyle='none'
 gg(2).Marker='p'
 gg(2).LineWidth = 0.8;
 gg(2).MarkerSize=1.5
xlabel(' Time Steps')
ax = gca;
ax.TickDir = 'out';
ylim([-26 5])
xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'r--','LineWidth',0.5,'HandleVisibility','off')
ax = gca; 
ax.TickDir = 'out';
title('Mode=3')
hold on


subplot(224)
gg=plot(0:length(a2_ww4_mat)-1,a2_ww4_mat',':')
 gg(1).Color = 'g'; 
 gg(5).Color = [0.9290 0.6940 0.1250];
 gg(3).Color = [0 0 1]; 
 gg(4).Color = [1 0 1];
 gg(2).Color = [0 0.4470 0.7410]; 
 gg(1).LineStyle='-'
 gg(5).LineStyle='--'
 gg(5).LineWidth = 0.8;
 gg(3).LineStyle=':'
 gg(3).LineWidth = 0.8;
 gg(4).LineStyle='-.'
gg(4).LineWidth = 0.8;
  gg(2).LineStyle='none'
 gg(2).Marker='p'
 gg(2).LineWidth = 0.8;
 gg(2).MarkerSize=1.5
xlabel(' Time Steps')
ax = gca; 
ax.TickDir = 'out';
ylim([-26 5])
xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'r--','LineWidth',0.5,'HandleVisibility','off')
ax = gca; 

ax.TickDir = 'out';
title('Mode=4')


figure(45)
subplot(211)

plot(0:length(a1_m_mat)-1, a1_m_mat', '-.o', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', '\rho_{1}^{1(k)}')
hold on
plot(0:length(a1_n_mat)-1, a1_n_mat', '--s', 'LineWidth', 1, 'MarkerSize', 3,'color','g', 'DisplayName', '\rho_{2}^{1(k)}')
plot(0:length(a1_p_mat)-1, a1_p_mat', ':p', 'LineWidth', 1, 'MarkerSize', 3,'color','[1 0 1]', 'DisplayName', '\rho_{3}^{1(k)}')
plot(0:length(a1_q_mat)-1, a1_q_mat', '-d', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', '\rho_{4}^{1(k)}')


legend('fontsize', 8)

xlabel('(a) Time Steps')
ylabel('Data-Based Case')
title('Area 1')


ax = gca;
ax.TickDir = 'out';


xlim([0, length(a1_q_mat)-1])
ylim([-1 22])
xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'k--','LineWidth',0.5,'HandleVisibility','off')


subplot(212)
plot(0:length(a2_m_mat)-1, a2_m_mat', '-.o', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', '\rho_{1}^{2(k)}')
hold on
plot(0:length(a2_n_mat)-1, a2_n_mat', '--s', 'LineWidth', 1, 'MarkerSize', 3, 'color','g','DisplayName', '\rho_{2}^{2(k)}')
plot(0:length(a2_p_mat)-1, a2_p_mat', ':p', 'LineWidth', 1, 'MarkerSize', 3,'color','[1 0 1]', 'DisplayName', '\rho_{3}^{2(k)}')
plot(0:length(a2_q_mat)-1, a2_q_mat', '-d', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', '\rho_{4}^{2(k)}')

legend('fontsize', 8)

xlabel('(b) Time Steps')
ylabel('Data-Based Case')
title('Area 2')

ax = gca;
ax.TickDir = 'out';

xlim([0, length(a2_q_mat)-1])
ylim([-1 30])
xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'k--','LineWidth',0.5,'HandleVisibility','off')



figure(46)
subplot(211)
plot(0:length(a1_p1_save)-1, a1_p1_save, '-.o', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 1');
hold on
plot(0:length(a1_p2_save)-1, a1_p2_save, '--s', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 2');
plot(0:length(a1_p3_save)-1, a1_p3_save, ':p', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 3');
plot(0:length(a1_p4_save)-1, a1_p4_save, '-d', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 4');

xlabel('(a)  Number of iterations')
ylabel('||P_{m}^{1(k)}-P_{m}^{1*}||')
title('Area 1')

legend('fontsize', 8)


ylim([-max([a1_p1_save, a1_p2_save, a1_p3_save, a1_p4_save])*0.2, max([a1_p1_save, a1_p2_save, a1_p3_save, a1_p4_save])*1.2])
xlim([0, 15])

ax = gca;
ax.TickDir = 'out';

xxxlim=get(gca,'Xlim');
hold on
plot(xxxlim,[0,0],'k--','LineWidth',0.5,'HandleVisibility','off')

subplot(212)
plot(0:length(a2_p1_save)-1, a2_p1_save, '-.o', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 1');
hold on
plot(0:length(a2_p2_save)-1, a2_p2_save, '--s', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 2');
plot(0:length(a2_p3_save)-1, a2_p3_save, ':p', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 3');
plot(0:length(a2_p4_save)-1, a2_p4_save, '-d', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 4');


xlabel('(b)  Number of iterations')
ylabel('||P_{m}^{2(k)}-P_{m}^{2*}||')
title('Area 2')

legend('fontsize', 8)

ylim([-max([a2_p1_save, a2_p2_save, a2_p3_save, a2_p4_save])*0.2, max([a2_p1_save, a2_p2_save, a2_p3_save, a2_p4_save])*1.2])
xlim([0, 15])


ax = gca;
ax.TickDir = 'out';

xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'k--','LineWidth',0.5,'HandleVisibility','off')



figure(47)
subplot(211)
plot(0:length(a1_k1_save)-1, a1_k1_save, '-.o', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 1');
hold on
plot(0:length(a1_k2_save)-1, a1_k2_save, '--s', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 2');
plot(0:length(a1_k3_save)-1, a1_k3_save, ':p', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 3');
plot(0:length(a1_k4_save)-1, a1_k4_save, '-d', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 4');


xlabel('(a)  Number of iterations')
ylabel('||K_{m}^{1(k)}-K_{m}^{1*}||')
title('Area 1')


legend('fontsize', 8)

ylim([-max([a1_k1_save, a1_k2_save, a1_k3_save, a1_k4_save])*0.2, max([a1_k1_save, a1_k2_save, a1_k3_save, a1_k4_save])*1.2])
xlim([0, 15])


ax = gca;
ax.TickDir = 'out';

xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'k--','LineWidth',0.5,'HandleVisibility','off')

subplot(212)
plot(0:length(a2_k1_save)-1, a2_k1_save, '-.o', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 1');
hold on
plot(0:length(a2_k2_save)-1, a2_k2_save, '--s', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 2');
plot(0:length(a2_k3_save)-1, a2_k3_save, ':p', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 3');
plot(0:length(a2_k4_save)-1, a2_k4_save, '-d', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 4');

xlabel('(b)  Number of iterations')
ylabel('||K_{m}^{2(k)}-K_{m}^{2*}||')
title('Area 2')


legend('fontsize', 8)


ylim([-max([a2_k1_save, a2_k2_save, a2_k3_save, a2_k4_save])*0.2, max([a2_k1_save, a2_k2_save, a2_k3_save, a2_k4_save])*1.2])
xlim([0, 15])

ax = gca;
ax.TickDir = 'out';

xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'k--','LineWidth',0.5,'HandleVisibility','off')


%% opt L

figure(48)
subplot(211)
plot(0:length(a1_l1_save)-1, a1_l1_save, '-.o', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 1');
hold on
plot(0:length(a1_l2_save)-1, a1_l2_save, '--s', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 2');
plot(0:length(a1_l3_save)-1, a1_l3_save, ':p', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 3');
plot(0:length(a1_l4_save)-1, a1_l4_save, '-d', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 4');

xlabel('(a)  Number of iterations')
ylabel('||L_{m}^{1(k)}-L_{m}^{1*}||')
title('Area 1')

legend('fontsize', 8)

ylim([-max([a1_l1_save, a1_l2_save, a1_l3_save, a1_l4_save])*0.2, max([a1_l1_save, a1_l2_save, a1_l3_save, a1_l4_save])*1.2])
xlim([0, 15])
ax = gca;
ax.TickDir = 'out';
xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'k--','LineWidth',0.5,'HandleVisibility','off')

subplot(212)
plot(0:length(a2_l1_save)-1, a2_l1_save, '-.o', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 1');
hold on
plot(0:length(a2_l2_save)-1, a2_l2_save, '--s', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 2');
plot(0:length(a2_l3_save)-1, a2_l3_save, ':p', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 3');
plot(0:length(a2_l4_save)-1, a2_l4_save, '-d', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Mode 4');

xlabel('(b)  Number of iterations')
ylabel('||L_{m}^{2(k)}-L_{m}^{2*}||')
title('Area 2')

legend('fontsize', 8)

ylim([-max([a2_l1_save, a2_l2_save, a2_l3_save, a2_l4_save])*0.2, max([a2_l1_save, a2_l2_save, a2_l3_save, a2_l4_save])*1.2])
xlim([0, 15])

ax = gca;
ax.TickDir = 'out';
xxxlim=get(gca,'Xlim'); 
hold on
plot(xxxlim,[0,0],'k--','LineWidth',0.5,'HandleVisibility','off')

%%   
figure(101)
legend_labels = {'x_1','x_2','x_3','x_4','x_5'};

subplot(221)
wg=plot(t1_save,x1_save,'Linewidth',1.5)
wg(1).LineStyle='--'
wg(1).Color='b'
wg(2).LineStyle='-'
wg(2).Color='r'
wg(3).LineStyle=':'
wg(3).Color='c'
wg(4).LineStyle='-.'
wg(4).Color='k'
wg(5).LineStyle='-'
wg(5).Color='[0 0.4470 0.7410]'

legend('x_1','x_2','x_3','x_4','x_5','location','Best','NumColumns',2,'EdgeColor', 'none')
 xlabel('(a)  Time (sec)')
 ylabel('State Response in Area 1')

ax = gca; 
ax.TickDir = 'out';
set(gca,'ytick') 
set(gca,'xtick') 


subplot(223)
wg=plot(t2_save,x2_save,'Linewidth',1.5)
wg(1).LineStyle='--'
wg(1).Color='b'
wg(2).LineStyle='-'
wg(2).Color='r'
wg(3).LineStyle=':'
wg(3).Color='c'
wg(4).LineStyle='-.'
wg(4).Color='k'
wg(5).LineStyle='-'
wg(5).Color='[0 0.4470 0.7410]'

legend('x_1','x_2','x_3','x_4','x_5','location','Best','NumColumns',2,'EdgeColor', 'none')
xlabel('(c)  Time (sec)')
ylabel('State Response in Area 2')

ax = gca; 

ax.TickDir = 'out';

set(gca,'ytick') 
set(gca,'xtick') 


subplot(222) 
wh=plot(t1_save,xa1_save,'Linewidth',1.5)
wh(1).LineStyle='--'
wh(1).Color='b'
wh(2).LineStyle='-'
wh(2).Color='r'
wh(3).LineStyle=':'
wh(3).Color='c'
wh(4).LineStyle='-.'
wh(4).Color='k'
wh(5).LineStyle='-'
wh(5).Color='[0 0.4470 0.7410]'
legend('x_1','x_2','x_3','x_4','x_5','location','Best','NumColumns',2,'EdgeColor', 'none')
xlabel('(b)  Time (sec)') 
ylabel('State Response in Area 1')
 ax = gca; 
ax.TickDir = 'out';
set(gca,'ytick') 
set(gca,'xtick') 
hold on

subplot(224) 
wh=plot(t2_save,xa2_save,'Linewidth',1.5)
wh(1).LineStyle='--'
wh(1).Color='b'
wh(2).LineStyle='-'
wh(2).Color='r'
wh(3).LineStyle=':'
wh(3).Color='c'
wh(4).LineStyle='-.'
wh(4).Color='k'
wh(5).LineStyle='-'
wh(5).Color='[0 0.4470 0.7410]'
legend('x_1','x_2','x_3','x_4','x_5','location','Best','NumColumns',2,'EdgeColor', 'none')
xlabel('(d)  Time (sec)') 
ylabel('State Response in Area 2')

 ax = gca; 
ax.TickDir = 'out';
set(gca,'ytick') 
set(gca,'xtick')

hold on






% 
figure(111)
subplot(221)
wh1=plot(t1_save,u1_save,'Linewidth',1.5)
wh1(1).LineStyle='--'
wh1(1).Color='r'
legend('controlled u','location','Best')
xlabel('(a) Time (sec)') 
ylabel('Control Input')
title('Area 1')

subplot(222)
wh2=stairs(t1_save,w1_save,'Linewidth',1.5)
wh2(1).LineStyle='--'
wh2(1).Color='b'
wh2(2).LineStyle='-'
wh2(2).Color='r'
legend('controlled w1','controlled w2','location','Best')
xlabel('(b) Time (sec)') 
ylabel('Disturbance Input')
title('Area 1')


subplot(223)
wh1=plot(t2_save,u2_save,'Linewidth',1.5)
wh1(1).LineStyle='--'
wh1(1).Color='r'

legend('controlled u','location','Best')
xlabel('(c) Time (sec)') 
ylabel('Control Input')
title('Area 2')

subplot(224)
wh2=stairs(t2_save,w2_save,'Linewidth',1.5)
wh2(1).LineStyle='--'
wh2(1).Color='b'
wh2(2).LineStyle='-'
wh2(2).Color='r'
legend('controlled w1','controlled w2','location','Best')
xlabel('(d) Time (sec)') 
ylabel('Disturbance Input')
title('Area 2')
hold on











