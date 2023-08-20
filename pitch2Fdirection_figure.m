% 09/02/2023
% generate figure pitch to force direction fitting, for wasp kinematics article
% stand along code
%% data input
Fx=[3.11600000000000,30.3060000000000,-29.5120000000000,-34.2930000000000,-40.0020000000000,1.34800000000000,8.05600000000000];
Fz=[57.3290000000000,60.1400000000000,80.4950000000000,92.5440000000000,75.5000000000000,62.5990000000000,60.1410000000000];
pitch=[37.6000000000000,29.0500000000000,65.2000000000000,67.3100000000000,66.1200000000000,35.6000000000000,36.6000000000000];
Vxy=[0.0910000000000000,0.0980000000000000,0.144000000000000,0.157000000000000,0.149000000000000,0.195000000000000,0.189000000000000];
%% plot data points
figure; hold on; axis equal; grid on;
relative_Vxy=Vxy-min(Vxy); %normalize body horizonta velocity so that minimal value is zero
relative_Vxy=relative_Vxy./max(relative_Vxy); %normalize body horizonta velocity so that maximal value is one
for cycle=1:length(Vxy)
    cycle_color=[relative_Vxy(cycle) 0 1-relative_Vxy(cycle)]; %set color to be blue at begging, red at ending
    scatter(pitch(cycle),atan2d(Fz(cycle),Fx(cycle)),50,cycle_color,'filled')
end
%% draw fit trendline
p = polyfit(pitch, atan2d(Fz,Fx), 1);
px = [min(pitch) max(pitch)];
py = polyval(p, px);
plot(px, py, 'LineWidth', 3,'color','k');
%% figure annotations
xlim([0 90])
ylim([0 180])
xlabel('Body pitch (degree)','fontsize',14)
ylabel('Force direction (degree)','fontsize',14)