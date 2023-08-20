%=========================================================================%
% 29/12/2021 By Amir Sarig                                                %
%                                                                         %
% - Plot 2D or 3D wing trajectory with body possition and chord direction.%
% - Plot is normilized to body length                                     %
% - Option of adding error bars for trajectory X & Z location (2D only)   %
%                                                                         %
% Called from 'wing_trajectory_data'                                      %
% Calls to func. 'wingtip_velocity_plot'                                        %
%=========================================================================%
function wing_trajectory_plot(tcorr,mean_x,mean_y,mean_z,V_tip_storage)
%%will be used in 2D plots
    [file,path] = uigetfile('*.*','Select insect profile image','D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Digitizing\Kinematics'); %ask user for side image of the flyer
    I = imread(fullfile(path,file)); %load image as variable
    existing_wingbase_XY_coord = questdlg('Use existing wingbase XY coord?','Existing wingbase XY coord?','Yes','No','No'); %ask user whether to record annimation
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
%% plot
    %ask user what to plot
%     plot_type=input('Plot 2D or 3D wingtip trajectory? (2/3)   ');
    body_length=sqrt((mean(mean(mean_x(1).data(:,:),1))-mean(mean(mean_x(2).data(:,:),1)))^2+(mean(mean(mean_y(1).data(:,:),1))-mean(mean(mean_y(2).data(:,:),1)))^2+(mean(mean(mean_z(1).data(:,:),1))-mean(mean(mean_z(2).data(:,:),1)))^2); %calc mean 3D body length
    Lbase_position=[mean(mean(mean_x(3).data(:,:),1)),mean(mean(mean_y(3).data(:,:),1)),mean(mean(mean_z(3).data(:,:),1))]; %mean 3D left wing base position
    Rbase_position=[mean(mean(mean_x(5).data(:,:),1)),mean(mean(mean_y(5).data(:,:),1)),mean(mean(mean_z(5).data(:,:),1))]; %mean 3D right wing base position
    for digitizing_point=1:size(mean_x,2) %repeat for each digiting point
        normalize(digitizing_point).coords=[(mean(mean_x(digitizing_point).data(:,:),1)-Lbase_position(1))/body_length; (mean(mean_y(digitizing_point).data(:,:),1)-Lbase_position(2))/body_length; (mean(mean_z(digitizing_point).data(:,:),1)-Lbase_position(3))/body_length]; %mean 3D: (digitizing_point minus wing-base) divide in body length. %mean is calculated along columns
    end
    for dimension=1:3 %calc mean of left&right wing. repeat for each axis
        mean_L_R_wingtip(dimension,:)=mean([normalize(4).coords(dimension,:); normalize(6).coords(dimension,:)],1); %mean is calculated along columns
        mean_L_R_chord_lead(dimension,:)=mean([normalize(7).coords(dimension,:); normalize(8).coords(dimension,:)],1); %mean is calculated along columns
        mean_L_R_chord_trail(dimension,:)=mean([normalize(9).coords(dimension,:); normalize(10).coords(dimension,:)],1); %mean is calculated along columns
        mean_L_R_chord_shift(dimension,:)=mean_L_R_wingtip(dimension,:)-mean([mean_L_R_chord_lead(dimension,:); mean_L_R_chord_trail(dimension,:)],1); %distance between middle of the cord to wingtip. mean is calculated along columns
    end
%     if plot_type==2
    figure
    hold on
        plot(mean_L_R_wingtip(1,:),mean_L_R_wingtip(3,:),'Color',[0.19855 0.7214 0.63095],'LineWidth',3) %plot mean 2D trajectory. mean is calculated along columns
%                     for ind=2:(size(coords.x_axis(pts2plot(i)).data,2)/5) %calc std of every 5th point, from the 2nd point
%                         errX(ind)=std(coords.x_axis(pts2plot(i)).data(:,ind));
%                         errZ(ind)=std(coords.z_axis(pts2plot(i)).data(:,ind));
%                     end %and draw it
    %                 h=errorbar(mean(mean_x(pts2plot(i)).data(:,2:5:end),1),mean(mean_z(pts2plot(i)).data(:,2:5:end),1),errX,errX,errZ,errZ,'.','Color',[0.19855 0.7214 0.63095]); %add error bars on the X & Y axis
                      %% Unsuccecful trial to draw shaded SD
    %                   fillX=[(mean(mean_x(pts2plot(i)).data(:,:),1)+errX),fliplr(mean(mean_x(pts2plot(i)).data(:,:),1)-errX)];
    %                   fillZ=[(mean(mean_z(pts2plot(i)).data(:,:),1)+errZ),fliplr(mean(mean_z(pts2plot(i)).data(:,:),1)-errZ)];
    %                   h=fill(fillX,fillZ,'r');
    %                   set(h,'alpha',.5);

    % ploting wing chord
        for chord=1:2:length(tcorr) %repeat for each 2nd chord data in the cycle
        %plot 
            plot([mean_L_R_chord_lead(1,chord), mean_L_R_chord_trail(1,chord)]+mean_L_R_chord_shift(1,chord),[mean_L_R_chord_lead(3,chord), mean_L_R_chord_trail(3,chord)]+mean_L_R_chord_shift(3,chord),'Color',[0.2081 0.1663 0.5292],'LineWidth',2) %plot mean 2D wing chord axis, shifted to wingtip location.
            plot(mean_L_R_chord_lead(1,chord)+mean_L_R_chord_shift(1,chord),mean_L_R_chord_lead(3,chord)+mean_L_R_chord_shift(3,chord),'o','Color',[0.2081 0.1663 0.5292],'MarkerFaceColor',[0.2081 0.1663 0.5292]) %plot mean 2D leading edge location as full circle, shifted to wingtip location. 
        end      
        h = image(insect_corners_x,-insect_corners_y,I); %insert image to the figure
        uistack(h,'bottom') %set flyer image as background
        xlabel('Longitudinal axis (body lengths)', 'FontSize',24)
        ylabel('Dorsoventral axis (body lengths)', 'FontSize',24)
        grid on
%     elseif plot_type==3
    set(gca,'FontSize',20)    
    figure
    hold on
            plot3(normalize(1).coords(1,:),normalize(1).coords(2,:),normalize(1).coords(3,:),'*','Color','k') %plot head mean 3D trajectory as astriks
            plot3([normalize(1).coords(1,:) normalize(2).coords(1,:)],[normalize(1).coords(2,:) normalize(2).coords(2,:)],[normalize(1).coords(3,:) normalize(2).coords(3,:)],'Color','k') %plot mean body longitude axis
            plot3(normalize(4).coords(1,:),normalize(4).coords(2,:),normalize(4).coords(3,:),'Color',[0.19855 0.7214 0.63095],'LineWidth',3) %plot left wingtip mean 3D trajectory
            plot3(normalize(6).coords(1,:),normalize(6).coords(2,:),normalize(6).coords(3,:),'Color',[0 1 0],'LineWidth',3) %plot right wingtip mean 3D trajectory
            for chord=1:2:length(tcorr) %repeat for each chord data in the cycle
                %left wing
                shift=normalize(4).coords(:,chord)-mean([normalize(7).coords(:,chord), normalize(9).coords(:,chord)],2); %caclculate distance between chord center to wing tip
                plot3([normalize(7).coords(1,chord) normalize(9).coords(1,chord)]+shift(1),[normalize(7).coords(2,chord) normalize(9).coords(2,chord)]+shift(2),[normalize(7).coords(3,chord) normalize(9).coords(3,chord)]+shift(3),'Color',[0.2081 0.1663 0.5292],'LineWidth',2) %plot mean 3D wing chord axis, shifted to wingtip location. mean is calculated along columns
                plot3(normalize(7).coords(1,chord)+shift(1), normalize(7).coords(2,chord)+shift(2), normalize(7).coords(3,chord)+shift(3),'o','Color',[0.2081 0.1663 0.5292],'MarkerFaceColor',[0.2081 0.1663 0.5292]) %plot mean 3D leading edge location as full circle, shifted to wingtip location. mean is calculated along columns
                clearvars shift
                %right wing
                shift=normalize(6).coords(:,chord)-mean([normalize(8).coords(:,chord), normalize(10).coords(:,chord)],2);
                plot3([normalize(8).coords(1,chord) normalize(10).coords(1,chord)]+shift(1),[normalize(8).coords(2,chord) normalize(10).coords(2,chord)]+shift(2),[normalize(8).coords(3,chord) normalize(10).coords(3,chord)]+shift(3),'Color',[0.5725 0.1412 0.1569],'LineWidth',2) %plot mean 3D wing chord axis, shifted to wingtip location. mean is calculated along columns
                plot3(normalize(8).coords(1,chord)+shift(1), normalize(8).coords(2,chord)+shift(2), normalize(8).coords(3,chord)+shift(3),'o','Color',[0.5725 0.1412 0.1569],'MarkerFaceColor',[0.5725 0.1412 0.1569]) %plot mean 3D leading edge location as full circle, shifted to wingtip location. mean is calculated along columns
                clearvars shift
            end
%             end         
        xlabel('Longitudinal axis (body lengths)', 'FontSize',24)
        ylabel('Left-Right axis (body lengths)', 'FontSize',24)
        zlabel('Dorsoventral axis (body lengths)', 'FontSize',24)
        xlim([-1.5 1.5])
        ylim([-2 2])
        view(100,30)
        flapping_movie(normalize) %create movie of flapping cycle
%     end
    set(gca,'FontSize',20)
    wingtip_velocity_plot(mean_L_R_wingtip,V_tip_storage,I,insect_corners_x,insect_corners_y)
    figure(2) %make wingtip trajectory figure active
end
%figure title is added in after the function is terminated