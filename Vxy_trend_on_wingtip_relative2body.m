%===========================================================================================%
%09.03.2022 By Amir Sarig
% - Plot each cycle's wingtip trajectory relative wing-base (color coded by body's mean Vxy during the cycle, cycles are ploted one over the other)
% - Plot cycles with max, median & min Vtip individually, Colored according to wingtip instentanious velocity
% - Calc wingtip 3D distance, normalized to body length
% - Plot stroke plane amlitude (color coded by body's mean Vxy during the cycle)
% - Plot velocity_trend_plots: Wingtip 3D speed relative to body & Body 3D speed (color coded by body's mean Vxy during the cycle)
% Called from 'extract_data_wasp_forward_flight'
%===========================================================================================%

function [I,insect_corners_x,insect_corners_y,wingtip_travel_dist]=Vxy_trend_on_wingtip_relative2body(Vxy,frq,rotated_coords,AllData_V_tip_storage,amplitude)
%% prepare data for plot
relative_Vxy=Vxy-min(Vxy); %normalize body horizonta velocity so that minimal value is zero
relative_Vxy=relative_Vxy./max(relative_Vxy); %normalize body horizonta velocity so that maximal value is one
[~,p]=sort(Vxy); %'p' describes the rearrangement of Vxy's elements along the sorted dimension such that Vxy(p) = sorted data
% [file,path] = uigetfile({'*.tif;*.jpg'},'Select insect profile image','D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Digitizing\Kinematics\WaspKinematics'); %ask user for side image of the flyer
path=('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Digitizing\Kinematics\WaspKinematics\');
file=('Wasp_body_profile_figure for wing_trajectory_plot.jpg');
I = imread(fullfile(path,file)); %load image as variable
existing_wingbase_XY_coord = questdlg(sprintf('Use existing wingbase XY coord?\n (if uncertain - choose the default)'),'Load data?','Yes','No','Yes'); %ask user whether to record annimation
if strcmp(existing_wingbase_XY_coord,'Yes') %if user choose to make use of saved matlab variables file 
    load(fullfile(path,'wingbase_XY_coord')) %load user selected variables to workspace
else %otherwise, user must have choosen a image
    figure
    image(I); %show insect image to user
    axis equal
    uiwait(msgbox('Click on the wing base location'));
    [x,y,button] = ginput(1); %recieve wing base coordinates from user
    currentFolder = pwd;
    close %close insect image
    insect_corners_x=[-x/size(I,2) 1+(-x/size(I,2))]; %calculate the lower and upper X values of the insect image location inside the plot, assuming the image length is 1
    insect_corners_y=[-y/size(I,1) (size(I,1)/size(I,2))-(-y/size(I,1))]; %calculate the lower and upper Y values of the insect image location inside the plot, assuming the image length is 1
    cd(path)
    save('wingbase_XY_coord','insect_corners_x','insect_corners_y'); %save wingbase_XY_coord in order to makes repeats quicker and identical
    cd(currentFolder)
end
%% calc max & min wingtip 3D velocity among all data
min_Vtip=10; max_Vtip=0; %set inital data
for wing_cycle=1:(length(Vxy))
    V_tip.x=mean([AllData_V_tip_storage(wing_cycle).body.left.x; AllData_V_tip_storage(wing_cycle).body.right.x],1);
    V_tip.y=mean([AllData_V_tip_storage(wing_cycle).body.left.y; AllData_V_tip_storage(wing_cycle).body.right.y],1);
    V_tip.z=mean([AllData_V_tip_storage(wing_cycle).body.left.z; AllData_V_tip_storage(wing_cycle).body.right.z],1);
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
%% calc normalize wingtip trajectory
% velocity_plot=figure; 
velocity_trend_plot=figure; hold on; title('Wingtip 3D speed relative to body'); xlabel('Flap cycle fraction'); ylabel('3D wingtip speed (m/s)');
body_velocity_trend_plot=figure; hold on; title('Body 3D speed'); xlabel('Flap cycle fraction'); ylabel('3D body speed (m/s)');
amplitude_plot=figure; hold on;
trajectory_plot=figure; hold on;
for wing_cycle=1:(length(Vxy)) %for each wing cycle
    body_length=sqrt((mean(rotated_coords(wing_cycle).data.x_axis(1).data)-mean(rotated_coords(wing_cycle).data.x_axis(2).data))^2+(mean(rotated_coords(wing_cycle).data.y_axis(1).data)-mean(rotated_coords(wing_cycle).data.y_axis(2).data))^2+(mean(rotated_coords(wing_cycle).data.z_axis(1).data)-mean(rotated_coords(wing_cycle).data.z_axis(2).data))^2); %calc mean 3D body length
    Lbase_position=[mean(rotated_coords(wing_cycle).data.x_axis(3).data),mean(rotated_coords(wing_cycle).data.y_axis(3).data),mean(rotated_coords(wing_cycle).data.z_axis(3).data)]; %mean 3D left wing base position
    Rbase_position=[mean(rotated_coords(wing_cycle).data.x_axis(5).data),mean(rotated_coords(wing_cycle).data.y_axis(5).data),mean(rotated_coords(wing_cycle).data.z_axis(5).data)]; %mean 3D right wing base position
    mean_L_R_wingbase=[mean([Lbase_position(1),Rbase_position(1)]),mean([Lbase_position(2),Rbase_position(2)]),mean([Lbase_position(3),Rbase_position(3)])];
    for digitizing_point=1:size(rotated_coords(wing_cycle).data.x_axis,2) %repeat for each digiting point: head;abd;Lbase;Ltip;Rbase;Rtip;Llead;Rlead;Ltrail;Rtrail
        normalize(digitizing_point).coords=[(rotated_coords(wing_cycle).data.x_axis(digitizing_point).data(:,:)-mean_L_R_wingbase(1))/body_length; (rotated_coords(wing_cycle).data.y_axis(digitizing_point).data(:,:)-mean_L_R_wingbase(2))/body_length; (rotated_coords(wing_cycle).data.z_axis(digitizing_point).data(:,:)-mean_L_R_wingbase(3))/body_length]; %(digitizing_point minus wing-base) divide in body length.
    end
%     %% calc wingtip 3D distance
%     for frame=1:(length(rotated_coords(wing_cycle).data.x_axis(4).data)-1) %repeat for each frame from digitizing data, exept the last
%         left_wingtip_travel_dist(frame)=sqrt((rotated_coords(wing_cycle).data.x_axis(4).data(:,frame+1)-rotated_coords(wing_cycle).data.x_axis(4).data(:,frame)).^2+(rotated_coords(wing_cycle).data.y_axis(4).data(:,frame+1)-rotated_coords(wing_cycle).data.y_axis(4).data(:,frame)).^2+(rotated_coords(wing_cycle).data.z_axis(4).data(:,frame+1)-rotated_coords(wing_cycle).data.z_axis(4).data(:,frame)).^2); %Calc 3D distance traveled between each to adjecent frames
%         right_wingtip_travel_dist(frame)=sqrt((rotated_coords(wing_cycle).data.x_axis(6).data(:,frame+1)-rotated_coords(wing_cycle).data.x_axis(6).data(:,frame)).^2+(rotated_coords(wing_cycle).data.y_axis(6).data(:,frame+1)-rotated_coords(wing_cycle).data.y_axis(6).data(:,frame)).^2+(rotated_coords(wing_cycle).data.z_axis(6).data(:,frame+1)-rotated_coords(wing_cycle).data.z_axis(6).data(:,frame)).^2); %Calc 3D distance traveled between each to adjecent frames
%     end
%     wingtip_travel_dist(wing_cycle)=mean([sum(left_wingtip_travel_dist),sum(right_wingtip_travel_dist)]); %sum distance between all frames, and calc mean between right&left wings
%     clear left_wingtip_travel_dist right_wingtip_travel_dist
    for dimension=1:3 %calc mean of left&right wing position. repeat for each axis
        mean_L_R_wingtip(dimension,:)=mean([normalize(4).coords(dimension,:); normalize(6).coords(dimension,:)],1); %mean is calculated along columns
    end
    mean_x=mean_L_R_wingtip(1,:); %extract x axis position of wingtip cycle
    mean_z=mean_L_R_wingtip(3,:); %extract z axis position of wingtip cycle
    %% calc wingtip 3D distance, normalized to body length
    for frame=1:(size(mean_L_R_wingtip,2)-1) %repeat for each frame from digitizing data, exept the last
        temp_wingtip_travel_dist(frame)=norm(mean_L_R_wingtip(:,frame+1)-mean_L_R_wingtip(:,frame)); %Calc 3D distance traveled between each to adjecent frames
    end
    wingtip_travel_dist(wing_cycle)=sum(temp_wingtip_travel_dist); %sum distance between all frames, and calc mean between right&left wings
    clear temp_wingtip_travel_dist
 %% plot wingtip trajectory relative 2 body data
    for i=1:1:length(rotated_coords(wing_cycle).data.x_axis(1).data)-1 %for each phase of wing cycle (exept the last)
        plot(mean_x(i:i+1),mean_z(i:i+1),'Color',[relative_Vxy(wing_cycle) 0 1-relative_Vxy(wing_cycle)],'LineWidth',3) %plot mean 2D trajectory betweem this phase to the next one. color it according to Vxy
    end
%     plot(mean_x(i+1:1),mean_z(i+1:1),'Color',[relative_Vxy(wing_cycle) 0 1-relative_Vxy(wing_cycle)],'LineWidth',4) %plot mean 2D trajectory betweem last phase to the first one. color it according to Vxy
%% plot SP amlitude
    if ~isnan(amplitude{1,1}) %if data exist
        figure(amplitude_plot)
        ventral_stroke=round((amplitude{wing_cycle,4}/100)*length(mean_x));
        dorsal_stroke_begging=round((amplitude{wing_cycle,5}/100)*length(mean_x));
        dorsal_stroke_x=mean(mean_x(end-dorsal_stroke_begging:end-1));
        dorsal_stroke_z=mean(mean_z(end-dorsal_stroke_begging:end-1));
        plot([mean_x(ventral_stroke),dorsal_stroke_x],[mean_z(ventral_stroke),dorsal_stroke_z],'Color',[relative_Vxy(wing_cycle) 0 1-relative_Vxy(wing_cycle)],'LineWidth',3)
    end
%% plot velocity_trend_plot
    figure(velocity_trend_plot)
    V_tip.x=mean([AllData_V_tip_storage(wing_cycle).body.left.x; AllData_V_tip_storage(wing_cycle).body.right.x],1); %L&R wings mean 
    V_tip.y=mean([AllData_V_tip_storage(wing_cycle).body.left.y; AllData_V_tip_storage(wing_cycle).body.right.y],1); %L&R wings mean
    V_tip.z=mean([AllData_V_tip_storage(wing_cycle).body.left.z; AllData_V_tip_storage(wing_cycle).body.right.z],1); %L&R wings mean
    Vtip_3d=((V_tip.x.^2)+(V_tip.y.^2)+(V_tip.z.^2)).^(1/2);%calc 3D velocity
    Vtip_range=[min_Vtip max_Vtip]; %store minimal and maximal wingtip velocity (of all cycle) 
    Vtip=Vtip_3d-min_Vtip; %normalize wingtip velocity so that minimal value (of all cycle) is zero
    Vtip=Vtip./max_Vtip; %normalize wingtip velocity so that maximal value (of all cycle) is one
    plot((1/length(Vtip)):(1/length(Vtip)):1,Vtip_3d,'Color',[relative_Vxy(wing_cycle) 0 1-relative_Vxy(wing_cycle)])
    figure(body_velocity_trend_plot)
    B_V_3D=((AllData_V_tip_storage(wing_cycle).body_center.x.^2)+(AllData_V_tip_storage(wing_cycle).body_center.y.^2)+(AllData_V_tip_storage(wing_cycle).body_center.z.^2)).^(1/2);%calc 3D velocity
    plot((1/length(Vtip)):(1/length(Vtip)):1,B_V_3D,'Color',[relative_Vxy(wing_cycle) 0 1-relative_Vxy(wing_cycle)])
    figure(trajectory_plot)
    %% plot cycles with max, median & min Vtip individually
    if p(wing_cycle)==1 || p(wing_cycle)==round(length(p)/2) || p(wing_cycle)==(length(p)) %for cycles with max, median & min Vxy
%         figure(velocity_plot)
        switch p(wing_cycle)
            case 1 %slowest
%                 subplot(3,1,1); 
                figure; hold on;
            case round(length(p)/2) %median
%                 subplot(3,1,2); 
                figure; hold on
            case length(p) %fastest
%                 subplot(3,1,3); 
                figure; hold on
        end
        % plot 
        for i=1:1:numel(Vtip)-1 %for each phase of wing cycle (exept the last)
            plot(mean_x(i:i+1),mean_z(i:i+1),'Color',[Vtip(i),sin((Vtip(i))*3.14),1-Vtip(i)],'LineWidth',2) %plot mean 2D trajectory betweem this phase to the next one. color it according to wingtip velocity
            hold on 
            plot(mean_x(i:i+1),mean_z(i:i+1),'Color',[Vtip(i),sin((Vtip(i))*3.14),1-Vtip(i)],'marker','o','MarkerFaceColor',[Vtip(i),sin((Vtip(i))*3.14),1-Vtip(i)]) %plot dot on this phase. color it according to wingtip velocity
        end
        axis equal
        ylim([-1 1.5])
        xlim([-1.1 1.1])
        xlabel('$$\hat{x}$$ (body length)','Interpreter','Latex', 'FontSize',24)%Longitudinal axis (body lengths)
        ylabel('$$\hat{z}$$ (body length)','Interpreter','Latex', 'FontSize',24)%Dorsoventral axis (body lengths)
        xlabel('x (body length)','FontSize',24)
        ylabel('z (body length)','FontSize',24)
        t=sprintf('Body''s horizontal speed: %1.2f m/s. \nWing frq.: %3.0f',Vxy(wing_cycle),frq(wing_cycle)); title(t)
        figure(trajectory_plot)
    end
end
%% plot annotations of main figure
h = image(insect_corners_x,-insect_corners_y,I); %insert image to the figure
uistack(h,'bottom') %set flyer image as background
cbh = colorbar ; %Create Colorbar
cbh.Ticks = linspace(0, 101, 11);  %Create & update 11 ticks from minimal wingtip velocity to maximal
c=[(0:0.01:1); zeros(1,101); (1:-0.01:0)]'; %create matrix of RGB colormap: Low-Blue, High-Red
colormap(c); %apply colormap that fit to the plot range
TickLabels = linspace(min(Vxy), max(Vxy), 11) ; %create labels
for i=1:11
    TickLabels_string(i,:) = num2str(TickLabels(i),'%1.3f'); %modify labels values to string
end
cbh.TickLabels=TickLabels_string; %update labels in figure
xlabel('$$\hat{x}$$ (body length)','Interpreter','Latex', 'FontSize',24)%Longitudinal axis (body lengths)
ylabel('$$\hat{z}$$ (body length)','Interpreter','Latex', 'FontSize',24)%Dorsoventral axis (body lengths)
xlabel('x (body length)','FontSize',24)
ylabel('z (body length)','FontSize',24)
axis equal
ylim([-1 1.5])
xlim([-1.1 1.1])
grid on
set(gca,'FontSize',20);
% title('Vxy trend on wingtip trajectory relative to body center');
%%
figure(amplitude_plot)
h = image(insect_corners_x,-insect_corners_y,I); %insert image to the figure
uistack(h,'bottom') %set flyer image as background
xlabel('$$\hat{x}$$ (body length)','Interpreter','Latex', 'FontSize',24)%Longitudinal axis (body lengths)
ylabel('$$\hat{y}$$ (body length)','Interpreter','Latex', 'FontSize',24)%Dorsoventral axis (body lengths)
xlabel('x (body length)','FontSize',24)
ylabel('z (body length)','FontSize',24)
xlim([-1.1 1.1])
% ylim([-1.5 1.5])
axis equal
grid on
set(gca,'FontSize',20);
%%
figure(velocity_trend_plot) %change active figure to wingtip 3D speed figure
legend(num2str(Vxy,2)); %add legend with Vxy values as repeat's name
figure(body_velocity_trend_plot) %change active figure to body 3D speed figure
legend(num2str(Vxy,2)); %add legend with Vxy values as repeat's name
ylim([0 2.5])
end