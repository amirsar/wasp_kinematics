%===========================================================================================%
%By Amir Sarig
% - Plot comparison between each cycle's upstroke to downstorke AoA & choed velcity 
% Called from 'extract_data_wasp_forward_flight'
%===========================================================================================%
%% extract data
load('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Publications\Wasp kinematics article\Data\phases_duration'); %manually prepared variable
for wing_cycle=1:(length(angles_storage)/2)
    downstroke=[];
    upstroke=[];
    AoA_data(wing_cycle).Vxy=angles_storage(2*wing_cycle-1).Vxy; %exctract flapping cycle Vxy
    Vxy(wing_cycle,1)=AoA_data(wing_cycle).Vxy; %will be used in the second loop
    Vxy(wing_cycle,2)=AoA_data(wing_cycle).Vxy; %will be used in the second loop
    AoA_data(wing_cycle).data(:,1)=angles_storage(2*wing_cycle-1).abs_AoA'; %exctract flapping cycle AoA
    tmpVchord=((AllData_coords_storage(wing_cycle).V_mid_chord.left+AllData_coords_storage(wing_cycle).V_mid_chord.right)./2)'; %R&L wings mean mid-chord speed in 3 axis
    Vchord(wing_cycle).data(:,1)=((tmpVchord(:,1).^2)+(tmpVchord(:,2).^2)+(tmpVchord(:,3).^2).^(0.5)); %3D mid-chord speed
    %% insert phases
    if cycle_phases{wing_cycle,1}==angles_storage(2*wing_cycle-1).cycle_name &&  cycle_phases{wing_cycle,1}==AllData_coords_storage(wing_cycle).cycle_name %verify that it is the same cycle
        [M,I]=max(mean([angles_storage(2*wing_cycle).elevation(1:10); angles_storage((2*wing_cycle)-1).elevation(1:10)],1)); %max elevation (mean of L&R wings) at the 1st half of the stroke
        cycle_phases{wing_cycle,9}=(I/length(angles_storage(2*wing_cycle).elevation))*100; %convert index to precentage
        phases=cycle_phases{wing_cycle,6:9}; %extract phases in percents
        phases_frames=cycle_phases{wing_cycle,3:5}; %extract phases in frames
        Vchord(wing_cycle).Vxy=AoA_data(wing_cycle).Vxy; %store Vxy in Vchord as well
    end
    wingtip_coords=(AllData_coords_storage(wing_cycle).original(4).data+AllData_coords_storage(wing_cycle).original(6).data)./2; %extract R&L wingtip mean coordinates
    amplitude(wing_cycle)=norm(wingtip_coords(length(wingtip_coords)-phases_frames(3)+1,:)-wingtip_coords(phases_frames(1),:)); %calc amplitude as clap minus ventral reversal point
    cycle_duration=48; %it's 50 frames, but include 1 frame (=2 extrapulated dat points) from the following cycle
%     [value,min_AoA_loc]=min(AoA_data(wing_cycle).data(1:round((phases(1)/100)*cycle_duration),1)); %find local min AoA during the downstroke
%     [value,max_AoA_loc]=max(AoA_data(wing_cycle).data(1:min_AoA_loc,1)); %find local min AoA during the downstroke
%     fling_phase(wing_cycle,1)=(max_AoA_loc/cycle_duration)*100; %store
    for frame=1:length(angles_storage(2*wing_cycle-1).abs_AoA)
        if frame<=(phases(4)/100)*cycle_duration
            AoA_data(wing_cycle).data(frame,2)=1; %'fling';
            Vchord(wing_cycle).data(frame,2)=1; %'fling';
        elseif frame<=(phases(1)/100)*cycle_duration
            AoA_data(wing_cycle).data(frame,2)=2; %'downstroke';
            Vchord(wing_cycle).data(frame,2)=2; %'downstroke';
        elseif frame<=((phases(1)+phases(2))/100)*cycle_duration
            AoA_data(wing_cycle).data(frame,2)=3; %'upstroke';
            Vchord(wing_cycle).data(frame,2)=3; %'upstroke';
        elseif frame<=cycle_duration
            AoA_data(wing_cycle).data(frame,2)=4; %'clap_break';
            Vchord(wing_cycle).data(frame,2)=4; %'clap_break';
        else
            AoA_data(wing_cycle).data(frame,2)=4; %'next_cycle';
            Vchord(wing_cycle).data(frame,2)=4; %'next_cycle';
        end
    end
    clear wingtip_coords
end
amplitude=table(amplitude', 'VariableNames', {'amplitude'});
amplitude=[cycle_phases(:,1:2) amplitude];
% fling_phase=table(round(fling_phase), 'VariableNames', {'fling_prcnt'});
% cycle_phases=[cycle_phases fling_phase];
% cycle_phases{:,6}=cycle_phases{:,6}-cycle_phases{:,9}; %correcting downstroke duration by substracting fling duration
%% calc
for i=1:2 %perform investigation for AoA & Vchord
    if i==1
        data=AoA_data;
    else
        data=Vchord;
    end
    figure
    downstroke_data=[];
    upstroke_data=[];
    for wing_cycle=1:(length(angles_storage)/2) %repeat for as many cycles in data
        downstroke=[];
        upstroke=[];
        [v,current_cycle]=min(Vxy(:,i)); %find cycle with min Vxy, so that data will extracted in increasing Vxy order
        Vxy(current_cycle,i)=1; %give large value so this cycle won't be choosen again
        downstroke=data(current_cycle).data(data(current_cycle).data(:,2)==2,1); %extract AoA from all downstroke frames in this cycle
        upstroke=data(current_cycle).data(data(current_cycle).data(:,2)==3,1); %extract AoA from all upstroke frames in this cycle
        %% t-test
        [h,p]=ttest2(downstroke,upstroke); %Two-sample t-test
        data(wing_cycle).down2upstroke_PV=p;
        subplot(1,(length(angles_storage)/2),wing_cycle)
        scatter(repmat(1,length(downstroke),1),downstroke)
        hold on
        scatter(repmat(2,length(upstroke),1),upstroke)
        xlim([0.5 2.5])
        s=sprintf('V=%1.3f, p=%1.2f',data(wing_cycle).Vxy,p);
        title(s,'FontSize',8)
        set(gca,'XTickLabel',{'downstroke','upstroke'})
        set(gca,'XTickLabelRotation',45)
        %% mean ratio
        data(wing_cycle).downstroke_mean=mean(downstroke);
        data(wing_cycle).upstroke_mean=mean(upstroke);
        data(wing_cycle).down2upstroke_MeanRation=mean(downstroke)/mean(upstroke);
        %% ANOVA data
        start_row=length(downstroke_data)+1;
        downstroke_data(start_row:(start_row-1+length(downstroke)),1)=downstroke;
        downstroke_data(start_row:(start_row-1+length(downstroke)),2)=repmat(data(current_cycle).Vxy,length(downstroke),1);
        start_row=length(upstroke_data)+1;
        upstroke_data(start_row:(start_row-1+length(upstroke)),1)=upstroke;
        upstroke_data(start_row:(start_row-1+length(upstroke)),2)=repmat(data(current_cycle).Vxy,length(upstroke),1);
    end
    if i==1
        AoA_data=data;
        AoA_downstroke_data=downstroke_data;
        AoA_upstroke_data=upstroke_data;
    else
        Vchord=data;
        Vchord_downstroke_data=downstroke_data;
        Vchord_upstroke_data=upstroke_data;
    end
    clear data 
end