% Draw scatter plots for the E.mundus kinematics article
clear
close all
load('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Publications\Wasp kinematics article\Data\reported_data.mat') %load kinematics data of each flapping cycle
load('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Publications\Wasp kinematics article\Data\wind_tunnel_mean_data.mat') %load wind tunnel data, mean of each individual
load('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Publications\Wasp kinematics article\Data\kinematics_mean_data.mat') %load kinematics data, mean of each individual
individuals=unique(report{:,1},'stable'); %extract serial number of the individual insects in the data 
for i=3:3 %each repeat create same set of scatter plots, but using different data: 1-. 2-. 3-latest version include only ventral facing ground flapping cycles kinematics, each cycle by itself.
    include_up_side_down=i-1;
    
    %% Pitch
    figure
    scatter(mean_windtunnel_data{:,5},mean_windtunnel_data{:,8},200,'bo','filled')
    marker=sprintf('pdsho');
    hold on
    switch include_up_side_down
        case 2
            for current_individual=1:length(individuals)
                flapping_cycles=individuals(current_individual)==report{:,1};
                scatter(report{flapping_cycles,7},report{flapping_cycles,17},200,'r',marker(current_individual),'LineWidth',3) 
                slope=58.3186; %value from statistica
                intercept=41.1791; %value from statistica
                h=refline(slope,intercept);
                set(h,'color','k')
            end
        case 1
            scatter(mean_kinematics_data{1:3,6},mean_kinematics_data{1:3,15},200,'ro','filled')
            scatter(mean_kinematics_data{4:5,6},mean_kinematics_data{4:5,15},200,'ro','LineWidth',2)
        case 0
            scatter(mean_kinematics_data{1:3,6},mean_kinematics_data{1:3,15},200,'ro','filled')
            slope=64.77654; %value from statistica
            intercept=38.83373; %value from statistica
            h=refline(slope,intercept);
            set(h,'color','k')
    end
    xlim([0 0.5])
    set(gca,'fontsize',24)
    ylabel('Body pitch ({\circ})','fontsize',35)
    ylim([0 90])
    legend({'Wind-tunnel';'Kinematics'})
    
    %% SP from body
    figure
    scatter(mean_windtunnel_data{:,5},mean_windtunnel_data{:,10},200,'bo','filled')
    hold on
    switch include_up_side_down
        case 2
            for current_individual=1:length(individuals)
                flapping_cycles=individuals(current_individual)==report{:,1};
                scatter(report{flapping_cycles,7},report{flapping_cycles,21},200,'r',marker(current_individual),'LineWidth',3)  
            end
        case 1
            scatter(mean_kinematics_data{1:3,6},mean_kinematics_data{1:3,20},200,'ro','filled')
            scatter(mean_kinematics_data{4:5,6},mean_kinematics_data{4:5,20},200,'ro','LineWidth',2)
        case 0
            scatter(mean_kinematics_data{1:3,6},mean_kinematics_data{1:3,20},200,'ro','filled')
    end
    xlim([0 0.5])
    set(gca,'fontsize',24)
    ylabel('Stroke plane ({\circ})','fontsize',35)
    ylim([0 135])
    legend({'Wind-tunnel';'Kinematics'})
    
    %% Pitch+SP
    figure
    % W_T_pitch_plus_SP=180-(mean_windtunnel_data{:,8}+mean_windtunnel_data{:,10});
    W_T_pitch_plus_SP=180-(mean_windtunnel_data{:,11});
    scatter(mean_windtunnel_data{:,5},W_T_pitch_plus_SP,200,'bo','filled')
    hold on
    % Kine_pitch_plus_SP=180-(mean_kinematics_data{:,15}+mean_kinematics_data{:,20});
    meanKine_pitch_plus_SP=180-(mean_kinematics_data{:,19});
    eachKine_pitch_plus_SP=180-(report{:,20});
    switch include_up_side_down
        case 2
            for current_individual=1:length(individuals)
                flapping_cycles=individuals(current_individual)==report{:,1};
                scatter(report{flapping_cycles,7},eachKine_pitch_plus_SP(flapping_cycles),200,'r',marker(current_individual),'LineWidth',3)  
                slope=-52.8002; %value from statistica
                intercept=37.8245; %value from statistica
                h=refline(slope,intercept);
                set(h,'color','k','linestyle','--') %nearly signifacnt
            end
        case 1
            scatter(mean_kinematics_data{1:3,6},meanKine_pitch_plus_SP(1:3),200,'ro','filled')
            scatter(mean_kinematics_data{4:5,6},meanKine_pitch_plus_SP(4:5),200,'ro','LineWidth',2)
        case 0
            scatter(mean_kinematics_data{1:3,6},meanKine_pitch_plus_SP(1:3),200,'ro','filled')
            slope=-46.2254; %value from statistica
            intercept=35.28466; %value from statistica
            h=refline(slope,intercept);
            set(h,'color','k','linestyle','--')
    end
    xlim([0 0.5])
    set(gca,'fontsize',24)
    ylabel('180-(Body pitch + SP)({\circ})','fontsize',35)
    % ylim([-10 90])
    legend({'Wind-tunnel';'Kinematics'})
    
    %% Frequency
    figure
    scatter(mean_windtunnel_data{:,5},mean_windtunnel_data{:,13},200,'bo','filled')
    hold on
    switch include_up_side_down
        case 2
            for current_individual=1:length(individuals)
                flapping_cycles=individuals(current_individual)==report{:,1};
                scatter(report{flapping_cycles,7},report{flapping_cycles,6},200,'r',marker(current_individual),'LineWidth',3)  
            end
        case 1
            scatter(mean_kinematics_data{1:3,6},mean_kinematics_data{1:3,5},200,'ro','filled')
            scatter(mean_kinematics_data{4:5,6},mean_kinematics_data{4:5,5},200,'ro','LineWidth',2)
        case 0
            scatter(mean_kinematics_data{1:3,6},mean_kinematics_data{1:3,5},200,'ro','filled')
    end
    xlim([0 0.5])
    set(gca,'fontsize',24)
    ylabel('Flapping frequency (Hz)','fontsize',35)
    ylim([0 350])
    legend({'Wind-tunnel';'Kinematics'})
    
    %% Amplitude distance
    figure
    scatter(mean_windtunnel_data{:,5},mean_windtunnel_data{:,14},200,'bo','filled')
    hold on
    switch include_up_side_down
        case 2
            for current_individual=1:length(individuals)
                flapping_cycles=individuals(current_individual)==report{:,1};
                scatter(report{flapping_cycles,7},report{flapping_cycles,27},200,'r',marker(current_individual),'LineWidth',3)
                slope=0.0007;
                intercept=0.0010;
                h=refline(slope,intercept);
                set(h,'color','k')
            end
        case 1
            scatter(mean_kinematics_data{1:3,6},mean_kinematics_data{1:3,24},200,'ro','filled')
            scatter(mean_kinematics_data{4:5,6},mean_kinematics_data{4:5,24},200,'ro','LineWidth',2)
        case 0
            scatter(mean_kinematics_data{1:3,6},mean_kinematics_data{1:3,24},200,'ro','filled')
            slope=0.000515;
            intercept=0.001031;
            h=refline(slope,intercept);
            set(h,'color','k')
    end
    xlim([0 0.5])
    set(gca,'fontsize',24)
    ylabel('Distance between stroke reversals (m)','fontsize',35)
    ylim([0.000 0.0015])
    legend({'Wind-tunnel';'Kinematics'})
end
figure %create figure with xlabel that will be attached to all of the above graphs
xlabel('Horizontal speed (m/s)','fontsize',35)