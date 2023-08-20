%Variable 'report' should be loaded to workspace
clearvars -except report
figure
hold on
color_palet=[0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4660, 0.6740, 0.1880; 0.4940, 0.1840, 0.5560; 0.6350, 0.0780, 0.1840];
color_counter=1;
for trial=1:size(report,1) %draw each trial's data
    if trial>1 %if it's not the 1st trial
        for i=1:size(color_palet,1) %check each color in the color palet
            if color_palet(i,4)==report{trial,1} %if it assign to the current individiual
                reuse_color=i; %use this color
                break %stop checking the rest of the colors
            else %if color wasn't assign with current individiual
                reuse_color=0; %flag as zero
            end
        end
        if reuse_color~=0; %if flag indicates that this individiual was assigned with color
            scatter(report{trial,6},report{trial,7},'markeredgecolor',color_palet(reuse_color,1:3),'MarkerFaceColor',color_palet(reuse_color,1:3)); %plot with assigned color
        else %if flag indicates that this individiual wasn't assigned with color
            color_palet(color_counter,4)=report{trial,1}; %assign new color with current individiual
            scatter(report{trial,6},report{trial,7},'markeredgecolor',color_palet(color_counter,1:3),'MarkerFaceColor',color_palet(color_counter,1:3)); %plot with assigned color
            color_counter=color_counter+1; %prepare new color
        end
    else %if it's the 1st trial
        color_palet(color_counter,4)=report{trial,1}; %assign new color with current individiual
        scatter(report{trial,6},report{trial,7},'markeredgecolor',color_palet(color_counter,1:3),'MarkerFaceColor',color_palet(color_counter,1:3)); %plot with assigned color
        color_counter=color_counter+1; %prepare new color
    end
    text(report{trial,6}+0.002,report{trial,7}+0.002,num2str(report{trial,2}),'Color','k') %add annotation of trial's name to the data
end
axis equal
xlim([0 0.3])
ylim([-0.009 0.04])
ylabel('Vertical Speed (m/s)')
xlabel('Horizontal Speed (m/s)')