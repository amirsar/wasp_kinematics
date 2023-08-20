function CFDbody4sim=CFDbody_measurments(wing_cycle,AllData_V_tip_storage_rawData,V_data,pitch,coords,CFDbody4sim)

CFDbody4sim(wing_cycle).cycle_name=AllData_V_tip_storage_rawData(wing_cycle).cycle_name; %store name
CFDbody4sim(wing_cycle).flapping_frq=V_data{1,4}; %store flapping frq.
CFDbody4sim(wing_cycle).velocity=[AllData_V_tip_storage_rawData(wing_cycle).body_center.x',AllData_V_tip_storage_rawData(wing_cycle).body_center.y',AllData_V_tip_storage_rawData(wing_cycle).body_center.z']; %store body center velocity vector
CFDbody4sim(wing_cycle).pitch=pitch.angle'; %store pitch angle vector
for point=1:length(coords(1).interpolation)
    CFDbody4sim(wing_cycle).body_length(point)=norm(coords(1).interpolation(:,point)-coords(2).interpolation(:,point)); %distance betweem head tip to abdomen tip
    CFDbody4sim(wing_cycle).body_width(point)=norm(coords(3).interpolation(:,point)-coords(5).interpolation(:,point)); %distance betweem wing bases
    CFDbody4sim(wing_cycle).Lwing_length(point)=norm(coords(3).interpolation(:,point)-coords(4).interpolation(:,point)); %distance betweem wing base tip wing tip
    CFDbody4sim(wing_cycle).Rwing_length(point)=norm(coords(5).interpolation(:,point)-coords(6).interpolation(:,point)); %distance betweem wing base tip wing tip
end
end