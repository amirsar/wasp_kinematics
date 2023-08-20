individuals=unique(report{:,1},'stable');
for current_individual=1:length(individuals)
    flapping_cycles=individuals(current_individual)==report{:,1};
    mean_kinematics_data(current_individual,:)=mean(report{flapping_cycles,:},1);
end
mean_kinematics_data = array2table(mean_kinematics_data,'VariableNames',report.Properties.VariableNames);
mean_kinematics_data(:,2)=[];
