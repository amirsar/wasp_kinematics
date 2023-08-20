function [clap_length,downstroke_length,upstroke_length] = cycle_phases_func(left_elev,right_elev,incidence)
    %% calc according to incidence
    clap_frames=find(incidence(:)>-10); %frames with elevation angle above -5 degrees, i.e. clap phase
    [reversal_point_value,reversal_point_frame]=min(incidence); %reversal_point is defined as the min. elevation
    clap_length=length(clap_frames)/length(incidence);
    downstroke_length=(reversal_point_frame-clap_frames(find(clap_frames<reversal_point_frame,1,'last')))/length(incidence); %(number of frames between reversal point to the last clap frame previous to the reversal point) devid to length of the cycle
    upstroke_length=(clap_frames(find(clap_frames>reversal_point_frame,1,'first'))-reversal_point_frame-1)/length(incidence); %(number of frames between reversal point to the first clap frame following to the reversal point, exluding the reversal point) devid to length of the cycle
    %% calc according to elevation
%     elevation=(left_elev+right_elev)./2;
%     clap=find(elevation(:)>75); %frames with elevation angle above 75 degrees, i.e. clap phase
%     [reversal_point_value,reversal_point_frame]=min(elevation); %reversal_point is defined as the min. elevation
%     clap_length=length(clap)/length(elevation);
%     if length(find(clap<reversal_point_frame, 1, 'last' ))>0 %if there is clap frames before reversal_point
%         downstroke_length=(reversal_point_frame-clap(find(clap<reversal_point_frame, 1, 'last' )))/length(elevation); %(number of frames between reversal point to the last clap frame previous to the reversal point) devid to length of the cycle
%     else
%         downstroke_length=(reversal_point_frame-0)/length(elevation); %(number of frames between reversal point begining of the cycle) devid to length of the cycle
%     end
%     upstroke_length=(clap(find(clap>reversal_point_frame, 1, 'first' ))-reversal_point_frame-1)/length(elevation); %(number of frames between reversal point to the first clap frame following to the reversal point, exluding the reversal point) devid to length of the cycle
end