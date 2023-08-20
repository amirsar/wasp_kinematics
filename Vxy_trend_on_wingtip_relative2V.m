%=========================================================================%
%15.03.2022 - By Amir Sarig
% - Plot each cycle's wingtip trajectory relative to body's velocity vector, color coded by body's mean Vxy during the cycle, cycles are ploted one over the other
% - Plot wingtip 3D speed relative to air
% - Plot mid-cord 3D velocity relative to air at body's velocity vector direction
% - Calculate pitch & store
% Called from 'extract_data_wasp_forward_flight'
% Calls to 'wingtip_world_viewpoint' which:
% - - Load & store insantious force data
% - - Plot each cycle's kinematics with force vectors overlay
%=========================================================================%
function [wingtip_travel_dist,Vtip_ground_3d,angles_storage,Force]=Vxy_trend_on_wingtip_relative2V(Vxy,coords,angles_storage,AllData_V_tip_storage)
%% calc max & min wingtip 3D velocity among all data
min_Vtip=10; max_Vtip=0; Force=[]; %set inital data
x_coords_limit=[0,0]; %will be used to set figure axis limits
for wing_cycle=1:(length(Vxy))
    V_tip.x=mean([AllData_V_tip_storage(wing_cycle).ground.left.x; AllData_V_tip_storage(wing_cycle).ground.right.x],1);
    V_tip.y=mean([AllData_V_tip_storage(wing_cycle).ground.left.y; AllData_V_tip_storage(wing_cycle).ground.right.y],1);
    V_tip.z=mean([AllData_V_tip_storage(wing_cycle).ground.left.z; AllData_V_tip_storage(wing_cycle).ground.right.z],1);
    Vtip=((median(V_tip.x,1).^2)+(median(V_tip.y,1).^2)+(median(V_tip.z,1).^2)).^(1/2);%calc 3D velocity
    if min(Vtip)<min_Vtip
        min_Vtip=min(Vtip);
    end
    if max(Vtip)>max_Vtip
        max_Vtip=max(Vtip);
    end
end
figure %create colorbar (in a different figure) for the wingtip velocity subplots
cbh = colorbar ; %Create Colorbar
cbh.Ticks = 0: 1/10: 1;  %Create & update 11 ticks from minimal wingtip velocity to maximal
c=[(0:0.01:1); sin((0:0.01:1)*3.14); (1:-0.01:0)]'; %create matrix of RGB colormap: Low-Blue, Middle-Green, High-Red
colormap(c); %apply colormap that fit to the plot range
TickLabels = linspace(min_Vtip, max_Vtip, 11) ; %create labels
for i=1:11
    TickLabels_string(i,:) = num2str(TickLabels(i),'%1.3f'); %modify labels values to string
end
cbh.TickLabels=TickLabels_string; %update labels in figure
%% prepare data for plot
relative_Vxy=Vxy-min(Vxy); %normalize body horizontal velocity so that minimal value is zero
relative_Vxy=relative_Vxy./max(relative_Vxy); %normalize body horizonta velocity so that maximal value is one
[~,p]=sort(Vxy); %'p' describes the rearrangement of Vxy's elements along the sorted dimension such that Vxy(p) = sorted data
r=0:1/(length(Vxy)-1):1; %'r' is a vector from 0 till 1 in the with the same amount of elements as Vxy
ranked_Vxy(p(:),1)=r; %'a' contains the rank of the corresponding element in Vxy
velocity_trend_plot=figure; hold on; title('Wingtip 3D speed relative to air'); xlabel('Flap cycle fraction'); ylabel('3D wingtip speed (m/s)');
% acceleration_trend_plot=figure; hold on; title('Wingtip 3D acceleration relative to air'); xlabel('Flap cycle fraction'); ylabel('3D wingtip acceleration (m/s^2)');
Wing_velocity_at_Body_velocity_vector_direction_trend_plot=figure; hold on; title('Mid-cord 3D velocity relative to air at body''s velocity vector direction'); xlabel('Flap cycle fraction'); ylabel('3D wingtip speed (m/s)');
tip_plot=figure; hold on
for wing_cycle=1:(length(Vxy)) %repeat for each wing cycle
    %% calc wingtip 3D distance
%     for frame=1:(length(coords(wing_cycle).interpolation(4).data)-1) %repeat for each frame from original digitizing data, exept the last
    for frame=1:(length(coords(wing_cycle).original(4).data)-1) %repeat for each frame from original digitizing data, exept the last
        left_wingtip_travel_dist(frame)=norm(coords(wing_cycle).original(4).data(frame+1,:)-coords(wing_cycle).original(4).data(frame,:)); %Calc 3D distance traveled between each to adjecent frames
        right_wingtip_travel_dist(frame)=norm(coords(wing_cycle).original(6).data(frame+1,:)-coords(wing_cycle).original(6).data(frame,:)); %Calc 3D distance traveled between each to adjecent frames
%         left_wingtip_travel_dist(frame)=norm(coords(wing_cycle).interpolation(4).data(:,frame+1)-coords(wing_cycle).interpolation(4).data(:,frame)); %Calc 3D distance traveled between each to adjecent frames
%         right_wingtip_travel_dist(frame)=norm(coords(wing_cycle).interpolation(6).data(:,frame+1)-coords(wing_cycle).interpolation(6).data(:,frame)); %Calc 3D distance traveled between each to adjecent frames
    end
    wingtip_travel_dist(wing_cycle)=mean([sum(left_wingtip_travel_dist),sum(right_wingtip_travel_dist)]); %sum distance between all frames, and calc mean between right&left wings
    clear left_wingtip_travel_dist right_wingtip_travel_dist
    %% calc normalize wingtip trajectory
    body_center=(coords(wing_cycle).interpolation(1).data(:,1)+coords(wing_cycle).interpolation(2).data(:,1))./2; %calc body center in the 1st frame
    for digitizing_point=1:size(coords(wing_cycle).interpolation,2) % start each cycle at (0,0,0), repeat for each digiting point. order: 1)head;2)abd;3)Lbase;4)Ltip;5)Rbase;6)Rtip;7)Llead;8)Rlead;9)Ltrail;10)Rtrail;
        normalize(digitizing_point).coords=(coords(wing_cycle).interpolation(digitizing_point).data-repmat(body_center,1,50))'; %(all frames of each digitizing_point minus body center at 1st frame) 
    end
%% calc rotated trajectory
    body_center=(normalize(1).coords(:,:)+normalize(2).coords(:,:))./2; %calc body center of each frame: (head+abd)/2
    velocity_vector=body_center(end,:)-body_center(1,:); %cycle's velocity_vector is BC's last frame minus 1st frame
    horizontal_angle=2*pi-atan2(velocity_vector(2),velocity_vector(1)); %calculate cycle's velocity_vector projection on XY plain
%     z_rotated=velocity_vector*[cos(horizontal_angle) -sin(horizontal_angle) 0; sin(horizontal_angle) cos(horizontal_angle) 0; 0 0 1];% rotated around vertical axis in order to have unified horizontal direction on the XZ plane. So that X is the horizontal direction of body velocity.
    z_rotated=([cos(horizontal_angle) -sin(horizontal_angle) 0; sin(horizontal_angle) cos(horizontal_angle) 0; 0 0 1]*velocity_vector')';% rotated around vertical axis in order to have unified horizontal direction on the XZ plane. So that X is the horizontal direction of body velocity.
    vertical_angle=2*pi-acos(dot(z_rotated([1 3]),[1,0])./(norm(z_rotated([1 3]))*norm([1,0]))); %calculate cycle's velocity_vector Vertical_angle from it's horizontal direction
%     V_y_rotated(wing_cycle,:)=z_rotated*[cos(-vertical_angle) 0 sin(-vertical_angle); 0 1 0; -sin(-vertical_angle) 0 cos(-vertical_angle)];% rotated in order to oriente XY plane
    V_y_rotated(wing_cycle,:)=([cos(-vertical_angle) 0 sin(-vertical_angle); 0 1 0; -sin(-vertical_angle) 0 cos(-vertical_angle)]*z_rotated')';% rotated in order to oriente XY plane
    for dim=1:3 %calc mean of left&right wing. repeat for each axis
        mean_L_R(1).coords(dim,:)=mean([normalize(4).coords(:,dim), normalize(6).coords(:,dim)],2); %wingtip mean is calculated along columns
        mean_L_R(2).coords(dim,:)=mean([normalize(7).coords(:,dim), normalize(8).coords(:,dim)],2); %leading edge mean is calculated along columns
        mean_L_R(3).coords(dim,:)=mean([normalize(9).coords(:,dim), normalize(10).coords(:,dim)],2); %trailing edge mean is calculated along columns
        mean_L_R(4).coords(dim,:)=normalize(1).coords(:,dim); %head
        mean_L_R(5).coords(dim,:)=normalize(2).coords(:,dim); %abdomen
    end
    mean_L_R(6).coords=((coords(wing_cycle).V_mid_chord.left'+coords(wing_cycle).V_mid_chord.right')./2)'; %mid-chord vector speed, L&R wings mean 
    for object=1:size(mean_L_R,2) %repeat for each object. order: 1)wingtip 2)lead 3)trail 4)head 5)abdomen 6)V_mid_chord
        for frame=1:size(mean_L_R(object).coords,2) %rotate each frame's wingtip coordinate by horizontal & vertical angles
            z_rotated=([cos(horizontal_angle) -sin(horizontal_angle) 0; sin(horizontal_angle) cos(horizontal_angle) 0; 0 0 1]*mean_L_R(object).coords(:,frame))';% rotated in order to oriente XZ plain
            y_rotated=([cos(-vertical_angle) 0 sin(-vertical_angle); 0 1 0; -sin(-vertical_angle) 0 cos(-vertical_angle)]*z_rotated')';% rotated in order to oriente XY plain
            object_relative2V(object).coords(:,frame)=y_rotated'; %store in a way that each axis is a row
        end
    end
    %% calc&store pitch relative2V & roll
    body_center=((object_relative2V(4).coords+object_relative2V(5).coords)./2); % rotated body center
    body_vector=(object_relative2V(4).coords-body_center)-(object_relative2V(5).coords-body_center); % normalized rotated head-abdomen
    wingtip=object_relative2V(1).coords-body_center; % normalized rotated wingtip
    clap_frame=round(0.9*size(object_relative2V(1).coords,2)); %middle of clap should be around 0.9 flapping cycle length
    transverse=cross(wingtip(:,clap_frame),body_vector(:,clap_frame)); %body transverse axis, from right to left of the body: cross(wings length,body_vector) 
    dorsal_norm=cross(body_vector(:,clap_frame),transverse); %dorsal perpendicular to the plane of the body
    if dorsal_norm(3)>=0 %if body's dorsal perpendicular is above body center
        pitch_sign=(1); %will leave pitch angle to range of 0 till 90. meaning ventral_face_ground
    else %if dorsal perpendicular is bellow body center
        pitch_sign=(-1); %will pass pitch angle to range of 90 till 180. meaning dorsal_face_ground
    end
    angles_storage(2*wing_cycle-1).pitch.angle2V=atan2d(body_vector(3,:),pitch_sign.*sqrt((body_vector(1,:).^2)+(body_vector(2,:).^2))); %calc body pitch relative to velocity vector. set the angles range according to body's dorsal direction at clap
    roll_angle=2*pi-acos(dot(dorsal_norm([1 3]),[0,1])./(norm(dorsal_norm([1 3]))*norm([0,1]))); %calculate cycle's roll_angle between Z axis to it's rotated body dorsal direction at clap
    for object=1:size(mean_L_R,2) %repeat for each object. order: 1)wingtip 2)lead 3)trail ...
        for frame=1:size(object_relative2V(object).coords,2) %rotate each frame's wingtip coordinate
            shifted=object_relative2V(object).coords(:,frame)-body_center(:,frame); %set body center as (0,0) because rotation is around the beginning of the axes
            x_rotated=([1 0 0; 0 cos(roll_angle) -sin(roll_angle); 0 sin(roll_angle) cos(roll_angle)]*shifted)';% rotated by roll_angle
            object_relative2V(object).coords(:,frame)=x_rotated'+body_center(:,frame); %store in a way that each axis is a row, return to original position
        end
    end
    if angles_storage(2*wing_cycle-1).pitch.ventral_face_down
        line_style='-'; %plot continous curve if ventral_face_down
    else
        line_style=':'; %plot dotted curve if ventral_face_up
    end
    for dim=1:3 %calc distance between wingtip to middle of chord. repeat for each axis
        mean_L_R_chord_shift(dim,:)=object_relative2V(1).coords(dim,:)-mean([object_relative2V(2).coords(dim,:); object_relative2V(3).coords(dim,:)],1); %distance between middle of the cord to wingtip. mean is calculated along columns
    end
    %% plot data
    % ploting wingtip
    figure(tip_plot)
    for i=1:1:size(object_relative2V(1).coords,2)-1 %for each phase of wing cycle (exept the last)
        plot(object_relative2V(1).coords(1,i:i+1),object_relative2V(1).coords(3,i:i+1),'Color',[relative_Vxy(wing_cycle) 0 1-relative_Vxy(wing_cycle)],'LineWidth',1.5,'LineStyle',line_style) %plot mean 2D trajectory betweem this phase to the next one. color it according to Vxy
    end
    if min(object_relative2V(1).coords(1,:))<x_coords_limit(1) %find minimal range for the figure
        x_coords_limit(1)=min(object_relative2V(1).coords(1,:)); %store it
    end
    if max(object_relative2V(1).coords(1,:))>x_coords_limit(2) %find maximal range for the figure
        x_coords_limit(2)=max(object_relative2V(1).coords(1,:)); %store it
    end   
    %% plot velocity_trend_plot
    figure(velocity_trend_plot)
    V_tip.x=mean([AllData_V_tip_storage(wing_cycle).ground.left.x; AllData_V_tip_storage(wing_cycle).ground.right.x],1); %L&R wings mean
    V_tip.y=mean([AllData_V_tip_storage(wing_cycle).ground.left.y; AllData_V_tip_storage(wing_cycle).ground.right.y],1); %L&R wings mean
    V_tip.z=mean([AllData_V_tip_storage(wing_cycle).ground.left.z; AllData_V_tip_storage(wing_cycle).ground.right.z],1); %L&R wings mean
    mean_x=object_relative2V(1).coords(1,:);
    mean_z=object_relative2V(1).coords(3,:);
    Vtip_ground_3d(wing_cycle,:)=((V_tip.x.^2)+(V_tip.y.^2)+(V_tip.z.^2)).^(1/2);%calc 3D velocity
    Vtip_range=[min_Vtip max_Vtip]; %store minimal and maximal wingtip velocity (of all cycle) 
    Vtip=Vtip_ground_3d(wing_cycle,:)-min_Vtip; %normalize wingtip velocity so that minimal value (of all cycle) is zero
    Vtip=Vtip./max_Vtip; %normalize wingtip velocity so that maximal value (of all cycle) is one
    plot((1/length(Vtip)):(1/length(Vtip)):1,Vtip_ground_3d(wing_cycle,:),'Color',[relative_Vxy(wing_cycle) 0 1-relative_Vxy(wing_cycle)])
    %% plot wingV@bodyVvector_direction_trend_plot
    figure(Wing_velocity_at_Body_velocity_vector_direction_trend_plot)
    plot((1/length(Vtip)):(1/length(Vtip)):1,object_relative2V(6).coords(1,:),'Color',[relative_Vxy(wing_cycle) 0 1-relative_Vxy(wing_cycle)])
    %% plot acceleration_trend_plot
%     figure(acceleration_trend_plot)
%     A_tip.x=mean([AllData_V_tip_storage(wing_cycle).acceleration.ground.left.x; AllData_V_tip_storage(wing_cycle).acceleration.ground.right.x],1);
%     A_tip.y=mean([AllData_V_tip_storage(wing_cycle).acceleration.ground.left.y; AllData_V_tip_storage(wing_cycle).acceleration.ground.right.y],1);
%     A_tip.z=mean([AllData_V_tip_storage(wing_cycle).acceleration.ground.left.z; AllData_V_tip_storage(wing_cycle).acceleration.ground.right.z],1);
%     Atip_ground_3d(wing_cycle,:)=((A_tip.x.^2)+(A_tip.y.^2)+(A_tip.z.^2)).^(1/2);%calc 3D velocity
%     plot((1/length(Vtip)):(1/length(Vtip)):1,Atip_ground_3d(wing_cycle,:),'Color',[relative_Vxy(wing_cycle) 0 1-relative_Vxy(wing_cycle)])
    %% plot cycles with max, median & min Vxy individually
%     if wing_cycle==2 || p(wing_cycle)==1 || p(wing_cycle)==round(length(p)/2) || p(wing_cycle)==(length(p)) %for cycles with max, median & min Vxy
%         switch p(wing_cycle)
%             case 1
%                 figure; hold on
%             case round(length(p)/2) 
%                 figure; hold on
%             case length(p)
%                 figure; hold on
%         end       
%         % plot wingtip
%         for i=1:1:numel(Vtip)-1 %for each phase of wing cycle (exept the last)
%             plot(mean_x(i:i+1),mean_z(i:i+1),'Color',[Vtip(i),sin((Vtip(i))*3.14),1-Vtip(i)],'LineWidth',2) %plot mean 2D trajectory betweem this phase to the next one. color it according to wingtip velocity
%             hold on
%             plot([object_relative2V(4).coords(1,i), object_relative2V(5).coords(1,i)],[object_relative2V(4).coords(3,i), object_relative2V(5).coords(3,i)],'Color','k','LineWidth',1) %plot body location.  
%             plot(object_relative2V(4).coords(1,i),object_relative2V(4).coords(3,i),'d','Color','k','MarkerFaceColor','k') %plot head location as full diamond.
%             scatter(mean_x(i:i+1),mean_z(i:i+1),[],[Vtip(i),sin((Vtip(i))*3.14),1-Vtip(i)],'filled') %plot dot on this phase. color it according to wingtip velocity
%         end       
%         % plot chord
%         for chord=1:2:size(object_relative2V(2).coords,2) %repeat for each 3nd chord data in the cycle
%             plot([object_relative2V(2).coords(1,chord), object_relative2V(3).coords(1,chord)]+mean_L_R_chord_shift(1,chord),[object_relative2V(2).coords(3,chord), object_relative2V(3).coords(3,chord)]+mean_L_R_chord_shift(3,chord),'Color',[0.2081 0.1663 0.5292],'LineWidth',1) %plot R&L wings mean 2D wing chord axis.
%             plot(object_relative2V(2).coords(1,chord)+mean_L_R_chord_shift(1,chord),object_relative2V(2).coords(3,chord)+mean_L_R_chord_shift(3,chord),'o','Color',[0.2081 0.1663 0.5292],'MarkerFaceColor',[0.2081 0.1663 0.5292]) %plot R&L wings mean 2D leading edge location as full circle, shifted to wingtip location. 
%         end
%         xlim([-0.0003 0.0018])
%         ylim([-0.0007 0.001])
%         axis equal
%         xlabel('$$\hat{x}$$ (m)','Interpreter','Latex', 'FontSize',24)%Body velocity axis
%         ylabel('$$\hat{y}$$ (m)','Interpreter','Latex', 'FontSize',24)%Vertical perpendicular axis
%         t=sprintf('body''s horizontal speed: %1.2f m/s',Vxy(wing_cycle)); title(t)
        Force=wingtip_world_viewpoint(wing_cycle,Vxy,coords,object_relative2V,Force); %plot new version
%     end
end
%% plot annotations of main figure
figure(tip_plot)
cbh = colorbar ; %Create Colorbar
cbh.Ticks = 0: 1/10: 1;  %Create & update 11 ticks from minimal wingtip velocity to maximal
c=[(0:0.01:1); zeros(1,101); (1:-0.01:0)]'; %create matrix of RGB colormap: Low-Blue, High-Red
colormap(c); %apply colormap that fit to the plot range
TickLabels = linspace(min(Vxy), max(Vxy), 11) ; %create labels
for i=1:11
    TickLabels_string(i,:) = num2str(TickLabels(i),'%1.3f'); %modify labels values to string
end
cbh.TickLabels=TickLabels_string; %update labels in figure
xlim(x_coords_limit.*1.2) %set x limit as 110% of minimal & maximal coords
axis equal
xlim([-0.0004 0.0016])
ylim([-9 11].*10^-4)
set(gca,'FontSize',20);
xlabel('$$\hat{x}$$ (m)','Interpreter','Latex', 'FontSize',24)%Body velocity axis
ylabel('$$\hat{z}$$ (m)','Interpreter','Latex', 'FontSize',24)%Vertical perpendicular axis
grid on
% title('Vxy trend on wingtip trajectory relative to body velocity');
figure(velocity_trend_plot) %change active figure to wingtip 3D speed figure
legend(num2str(Vxy,2)); %add legend with Vxy values as repeat's name
end