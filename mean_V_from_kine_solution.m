%===========================================================================================%
%21.11.2021
% - Calc&store for each wing cycle: 'Trial', 'Cycle', 'Advance ration', 'Frq.', 'Vxyz', 'Vxy', 'Vz','V_Ltip','V_Rtip' & interpolated coords of:head,abd,Ltip,Rtip
% called from 'AnglesExcel_Amir' & 'extract_data_wasp_forward_flight'
% calls: 'body_pitch' & 'velocity' functions
%===========================================================================================%
function [V_data,V_tip_storage,pitch,coords_storage] = mean_V_from_kine_solution(coords, data, continuous, start_frame,end_frame,Hz,V_data,ind,idx,draw_plot,tcorr,V_tip_storage,prior_frames,posterior_frames,elevation,pitch,coords_storage)
%% add frames in edges in order to avoid miscalculation using velocity calculation. They are revomed afterwards
if ~lower(continuous)=='f' %if cycles are continuous
    end_frame=end_frame-1; %because the givven cycles are overlapping
end
%% Obtaining body & wing points
head = coords(1).extended';
abd = coords(2).extended';
Lbase = coords(3).extended';
Ltip = coords(4).extended';
Rbase = coords(5).extended';
Rtip = coords(6).extended';
Llead = coords(7).extended';
Rlead = coords(8).extended';
Ltrail = coords(9).extended';
Rtrail = coords(10).extended';
Lmid_chord=(Llead'+Ltrail')./2; %calc mid chord at ~0.7 wing length
Rmid_chord=(Rlead'+Rtrail')./2; %calc mid chord at ~0.7 wing length
%% calc velocity. *because velocity func. enter error in edges, and because wingtip data is almost circular - adding to data 5 frames in each edge from the other edge
xyz=((head'+abd')./2); %of body center
Ltip=Ltip';
Rtip=Rtip';
v_body=velocity(xyz,Hz);
V_G_Ltip=velocity(Ltip,Hz);
V_G_Rtip=velocity(Rtip,Hz);
V_Lmid_chord=velocity(Lmid_chord,Hz);
V_Rmid_chord=velocity(Rmid_chord,Hz);
% V_G_Ltip4A=V_G_Ltip;
% V_G_Rtip4A=V_G_Rtip;
%removing data frames that were added to evoid velocity func. from entering errors
v_body=v_body((1+prior_frames):(end-posterior_frames),:);
V_G_Ltip=V_G_Ltip((1+prior_frames):(end-posterior_frames),:);
V_G_Rtip=V_G_Rtip((1+prior_frames):(end-posterior_frames),:);
V_B_Ltip=[V_G_Ltip(:,1)-v_body(:,1),V_G_Ltip(:,2)-v_body(:,2),V_G_Ltip(:,3)-v_body(:,3)];
V_B_Rtip=[V_G_Rtip(:,1)-v_body(:,1),V_G_Rtip(:,2)-v_body(:,2),V_G_Rtip(:,3)-v_body(:,3)];
V_Lmid_chord=V_Lmid_chord((1+prior_frames):(end-posterior_frames),:);
V_Rmid_chord=V_Rmid_chord((1+prior_frames):(end-posterior_frames),:);
%% interpollate cycle's wing tip & mid_chord speed to uniform data size
corrected_data_V_body=interpolate_data(v_body,tcorr);
corrected_data_V_G_Ltip=interpolate_data(V_G_Ltip,tcorr);
corrected_data_V_G_Rtip=interpolate_data(V_G_Rtip,tcorr);
corrected_data_V_B_Ltip=interpolate_data(V_B_Ltip,tcorr);
corrected_data_V_B_Rtip=interpolate_data(V_B_Rtip,tcorr);
corrected_data_V_Lmid_chord=interpolate_data(V_Lmid_chord,tcorr);
corrected_data_V_Rmid_chord=interpolate_data(V_Rmid_chord,tcorr);
%% acceleration
% corrected_data_V_G_Ltip4A=interpolate_data(V_G_Ltip4A,tcorr);
% corrected_data_V_G_Rtip4A=interpolate_data(V_G_Rtip4A,tcorr);
% s=1+prior_frames/length(xyz);
% e=(length(xyz)+posterior_frames)/length(xyz);
%% store data
row=size(V_tip_storage.body_center.x,1);
V_tip_storage.body_center.x(row+1,:)=corrected_data_V_body(1,:); %store X axis velocity of body center
V_tip_storage.body_center.y(row+1,:)=corrected_data_V_body(2,:); %store Y axis velocity of body center
V_tip_storage.body_center.z(row+1,:)=corrected_data_V_body(3,:); %store Z axis velocity of body center
V_tip_storage.ground.left.x(row+1,:)=corrected_data_V_G_Ltip(1,:); %store X axis velocity of left wing tip relative to air
V_tip_storage.ground.left.y(row+1,:)=corrected_data_V_G_Ltip(2,:); %store Y axis velocity of left wing tip relative to air
V_tip_storage.ground.left.z(row+1,:)=corrected_data_V_G_Ltip(3,:); %store Z axis velocity of left wing tip relative to air
V_tip_storage.ground.right.x(row+1,:)=corrected_data_V_G_Rtip(1,:); %store X axis velocity of right wing tip relative to air
V_tip_storage.ground.right.y(row+1,:)=corrected_data_V_G_Rtip(2,:); %store Y axis velocity of right wing tip relative to air
V_tip_storage.ground.right.z(row+1,:)=corrected_data_V_G_Rtip(3,:); %store Z axis velocity of right wing tip relative to air
% V_tip_storage.acceleration.ground.left.x(row+1,:)=corrected_data_A_G_Ltip(1,:); %store X axis acceleration of left wing tip relative to air
% V_tip_storage.acceleration.ground.left.y(row+1,:)=corrected_data_A_G_Ltip(2,:); %store Y axis acceleration of left wing tip relative to air
% V_tip_storage.acceleration.ground.left.z(row+1,:)=corrected_data_A_G_Ltip(3,:); %store Z axis acceleration of left wing tip relative to air
% V_tip_storage.acceleration.ground.right.x(row+1,:)=corrected_data_A_G_Rtip(1,:); %store X axis acceleration of right wing tip relative to air
% V_tip_storage.acceleration.ground.right.y(row+1,:)=corrected_data_A_G_Rtip(2,:); %store Y axis acceleration of right wing tip relative to air
% V_tip_storage.acceleration.ground.right.z(row+1,:)=corrected_data_A_G_Rtip(3,:); %store Z axis acceleration of right wing tip relative to air
V_tip_storage.body.left.x(row+1,:)=corrected_data_V_B_Ltip(1,:); %store X axis velocity of left wing tip relative to body
V_tip_storage.body.left.y(row+1,:)=corrected_data_V_B_Ltip(2,:); %store Y axis velocity of left wing tip relative to body
V_tip_storage.body.left.z(row+1,:)=corrected_data_V_B_Ltip(3,:); %store Z axis velocity of left wing tip relative to body
V_tip_storage.body.right.x(row+1,:)=corrected_data_V_B_Rtip(1,:); %store X axis velocity of right wing tip relative to body
V_tip_storage.body.right.y(row+1,:)=corrected_data_V_B_Rtip(2,:); %store Y axis velocity of right wing tip relative to body
V_tip_storage.body.right.z(row+1,:)=corrected_data_V_B_Rtip(3,:); %store Z axis velocity of right wing tip relative to body
row=size(coords_storage,1);
for dig_point=1:size(coords,2)
    coords_storage(row+1).interpolation(dig_point).data=coords(dig_point).interpolation;
    coords_storage(row+1).original(dig_point).data=coords(dig_point).original;
    coords_storage(row+1).extended(dig_point).data=coords(dig_point).extended;
    coords_storage(row+1).extended(dig_point).data=coords(dig_point).defualt;
end
coords_storage(row+1).prior_frames=prior_frames;
coords_storage(row+1).posterior_frames=posterior_frames;
coords_storage(row+1).V_mid_chord.left=corrected_data_V_Lmid_chord; %store velocity of left mid chord at 0.7 wing tip relative to air
coords_storage(row+1).V_mid_chord.right=corrected_data_V_Rmid_chord; %store velocity of right mid chord at 0.7 wing tip relative to air
%% calc mean speed
frq=1/((end_frame-start_frame)/Hz);
Vx=mean(v_body(:,1));
Vy=mean(v_body(:,2));
Vz=mean(v_body(:,3));
Vxy=mean(sqrt((v_body(:,1).^2)+(v_body(:,2).^2)));
Vxyz=mean(sqrt((v_body(:,1).^2)+(v_body(:,2).^2)+(v_body(:,3).^2)));
VLtip_ground=mean(sqrt((V_G_Ltip(:,1).^2)+(V_G_Ltip(:,2).^2)+(V_G_Ltip(:,3).^2)));
VRtip_ground=mean(sqrt((V_G_Rtip(:,1).^2)+(V_G_Rtip(:,2).^2)+(V_G_Rtip(:,3).^2)));
VLtip_body=mean(sqrt(((V_G_Ltip(:,1)-v_body(:,1)).^2)+((V_G_Ltip(:,2)-v_body(:,2)).^2)+((V_G_Ltip(:,3)-v_body(:,3)).^2)));
VRtip_body=mean(sqrt(((V_G_Rtip(:,1)-v_body(:,1)).^2)+((V_G_Rtip(:,2)-v_body(:,2)).^2)+((V_G_Rtip(:,3)-v_body(:,3)).^2)));
%% store data
[row,c]=size(V_data);
V_data{row+1,:}=[ind,idx,(Vxyz/mean([VLtip_body,VRtip_body])),frq,Vxyz,Vxy,Vx,Vy,Vz,VLtip_ground,VRtip_ground,VLtip_body,VRtip_body];
%*hovering=Advance ration<0.1 (Ellington 1985 part3 p.55)
if c==0 %define columns names
    V_data.Properties.VariableNames={'Trial', 'Cycle', 'Advance_ration', 'Frq', 'Vxyz', 'Vxy', 'Vx', 'Vy', 'Vz','V_Ltip_ground','V_Rtip_ground','V_Ltip_body','V_Rtip_body'};
end
pitch=body_pitch(pitch,ind,idx,row+1,head(:,(1+prior_frames):(end-posterior_frames)),abd(:,(1+prior_frames):(end-posterior_frames)),Ltip((1+prior_frames):(end-posterior_frames),:)',Rtip((1+prior_frames):(end-posterior_frames),:)',elevation,tcorr,draw_plot); %calc body pitch
if draw_plot
    plot_Vxy_Vs_Vz(V_data{:,1},V_data{:,2},V_data{:,6},V_data{:,9},data.config); %plot Horizontal Vs. Vertical Flight Speed
    grid on
    axis square
end
end