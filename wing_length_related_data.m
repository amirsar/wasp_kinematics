%===========================================================================================%
%By Amir Sarig
% -  Calculate wingtip's traveled distance normalized to body length
% Called from 'extract_data_wasp_forward_flight'
%===========================================================================================%
function distance=wing_length_related_data(AllData_coords_storage)
wasps_wing_length=[];
for wing_cycle=1:length(AllData_coords_storage) %repeaet for each wing cycle
    coords=AllData_coords_storage(wing_cycle).original; %extract origial digizing coords
    distance(wing_cycle).cycle_name=AllData_coords_storage(wing_cycle).cycle_name; %extract cycle serial name
    wing_base=(coords(3).data+coords(5).data)./2; %mean of R & L wing base loctaion
    wing_tip=(coords(4).data+coords(6).data)./2; %mean of R & L wing tip loctaion
    body_center=(coords(1).data+coords(2).data)./2; %mean between head to abdomen tip
    if isempty(wasps_wing_length) %if no wing length data were stored
        wasps=0; %insert zero in order to allow the following query
    else
        wasps=wasps_wing_length(:,1); %extract wasps number
    end
    if any(wasps==floor(distance(wing_cycle).cycle_name)) %if wing_length already calculated for the current wasp
        mean_wing_length=wasps_wing_length(wasps_wing_length(:,1)==floor(distance(wing_cycle).cycle_name),2); %exctract wing length of current wasp
    else %if wing_length for the current wasp wasn't found
        for frame=1:length(wing_tip) %repeaet for each frame
            wing_length(frame)=norm(wing_tip(frame,:)-wing_base(frame,:)); %calc size of wing length
        end
        mean_wing_length=mean(wing_length); %mean size in this cycle
        next_cycle=wing_cycle+1; %will be used in the following 'while' loop
        while floor(distance(wing_cycle).cycle_name)==floor(AllData_coords_storage(next_cycle).cycle_name) %repeat as long as the current cycle and the compared cycle belong to the same wasp
            tmp_wing_base=(AllData_coords_storage(next_cycle).original(3).data+AllData_coords_storage(next_cycle).original(5).data)./2; %mean of R & L wing base loctaion
            tmp_wing_tip=(AllData_coords_storage(next_cycle).original(4).data+AllData_coords_storage(next_cycle).original(6).data)./2; %mean of R & L wing tip loctaion
            for frame=1:length(tmp_wing_tip) %repeaet for each frame
                tmp_wing_length(frame)=norm(tmp_wing_tip(frame,:)-tmp_wing_base(frame,:)); %size of wing length
            end
            tmp_mean_wing_length(next_cycle-wing_cycle)=mean(tmp_wing_length); %mean size in this cycle. If the variable store data from several cycles, it store all
            next_cycle=next_cycle+1; %check the next cycle
            if next_cycle>length(AllData_coords_storage) %if next cycle doesn't exist
                break %terminate loop
            end
        end
        if wing_cycle-next_cycle>1 %if additional mean_wing_length were calculated for the current wasp 
            mean_wing_length=mean([mean_wing_length tmp_mean_wing_length]); %mean wing length from all cycles of the same wasp
        end
        wasps_wing_length(size(wasps_wing_length,1)+1,:)=[floor(distance(wing_cycle).cycle_name),mean_wing_length]; %store mean wing length of current wasp
    end
    normalized_wing_tip=(wing_tip-body_center)/mean_wing_length; %normalize wingtip location relative to wing base, and it's ratio to wing length
    for frame=2:length(normalized_wing_tip) %repeaet for each frame
        temp_distance(frame)=norm(normalized_wing_tip(frame,:)-normalized_wing_tip(frame-1,:)); %calc distacnce between adjecnt frames
    end
    distance(wing_cycle).distance=sum(temp_distance); %wingtip normalized traveled distance during wing cycle
    clearvars -except distance wasps_wing_length AllData_coords_storage
end
end