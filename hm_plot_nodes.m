function hm_plot_nodes
%
%
%

%% define areas
dpCo{1}=[300 500]*1e-9; %m
tmCo{1}=[0 11.5]; %hours

dpCo{2}=[300 500]*1e-9; %m
tmCo{2}=[11.5 24]; %hours

dpCo{3}=[100 300]*1e-9; %m
tmCo{3}=[0 7]; %hours

dpCo{4}=[100 300]*1e-9; %m
tmCo{4}=[7 15]; %hours

dpCo{5}=[100 300]*1e-9; %m
tmCo{5}=[15 24]; %hours

dpCo{6}=[20 100]*1e-9; %m
tmCo{6}=[0 7]; %hours

dpCo{7}=[20 100]*1e-9; %m
tmCo{7}=[7 13]; %hours

dpCo{8}=[50 100]*1e-9; %m
tmCo{8}=[13 24]; %hours

dpCo{9}=[30 50]*1e-9; %m
tmCo{9}=[13 16]; %hours

dpCo{10}=[30 50]*1e-9; %m
tmCo{10}=[16 20]; %hours

dpCo{11}=[30 50]*1e-9; %m
tmCo{11}=[20 24]; %hours

dpCo{12}=[3 20]*1e-9; %m
tmCo{12}=[0 5.5]; %hours

dpCo{13}=[3 20]*1e-9; %m
tmCo{13}=[5.5 10.5]; %hours

dpCo{14}=[3 20]*1e-9; %m
tmCo{14}=[10.5 13]; %hours

dpCo{15}=[10 30]*1e-9; %m
tmCo{15}=[13 19]; %hours

dpCo{16}=[10 30]*1e-9; %m
tmCo{16}=[19 24]; %hours

dpCo{17}=[3 10]*1e-9; %m
tmCo{17}=[13 17]; %hours

dpCo{18}=[3 10]*1e-9; %m
tmCo{18}=[17 24]; %hours

%% plot areas
xl=xlim;
hold on
dayStart=floor(mean(xl));
for i=1:length(dpCo)
    plot([dayStart+tmCo{i}(1)/24 dayStart+tmCo{i}(2)/24],[dpCo{i}(1) dpCo{i}(1)],'k')
    plot([dayStart+tmCo{i}(1)/24 dayStart+tmCo{i}(2)/24],[dpCo{i}(2) dpCo{i}(2)],'k')
    plot([dayStart+tmCo{i}(1)/24 dayStart+tmCo{i}(1)/24],[dpCo{i}(1) dpCo{i}(2)],'k')
    plot([dayStart+tmCo{i}(2)/24 dayStart+tmCo{i}(2)/24],[dpCo{i}(1) dpCo{i}(2)],'k')
end