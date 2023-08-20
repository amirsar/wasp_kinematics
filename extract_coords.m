function [data, coords] = extract_coords(data,start_frame,end_frame,prior_frames,posterior_frames,switch_wings,Hz,tcorr,use_filter)
%% extracting data
% Obtaining body & wing points, correct digitizing errors, filter data
% function is called from 'extract_coords' & 'AnglesExcel_Amir' scripts
% calling 'filter_digitizing_points' function
kine_frames = (start_frame-prior_frames):(end_frame+posterior_frames); %cycle frame range with addition of frames used for derevative calculation 
head = data.kine.body.data.coords(1,:,kine_frames); head = reshape(head,3,size(head,3)); %reshape so each column is a single dot
abd = data.kine.body.data.coords(2,:,kine_frames); abd = reshape(abd,3,size(abd,3)); 
Lbase = data.kine.left_wing.data.coords(1,:,kine_frames); Lbase = reshape(Lbase,3,size(Lbase,3)); 
Ltip = data.kine.left_wing.data.coords(2,:,kine_frames); Ltip = reshape(Ltip,3,size(Ltip,3)); 
Rbase = data.kine.right_wing.data.coords(1,:,kine_frames); Rbase = reshape(Rbase,3,size(Rbase,3));
Rtip = data.kine.right_wing.data.coords(2,:,kine_frames); Rtip = reshape(Rtip,3,size(Rtip,3)); 
Llead = data.kine.left_wing.data.model_coords(:,13,kine_frames); Llead = reshape(Llead,3,size(Llead,3)); %third point on the wing's leading edge
Rlead = data.kine.right_wing.data.model_coords(:,13,kine_frames); Rlead = reshape(Rlead,3,size(Rlead,3)); 
Ltrail = data.kine.left_wing.data.model_coords(:,31,kine_frames); Ltrail = reshape(Ltrail,3,size(Ltrail,3)); %third point on the wing's trailing edge
Rtrail = data.kine.right_wing.data.model_coords(:,31,kine_frames); Rtrail = reshape(Rtrail,3,size(Rtrail,3)); 

%% Digitizing corretion
switch switch_wings(1) %correction manipulation on digitizion
    case 1 %Digitizing is fine
        %no need in more actions
    case 2 %wings should be switched
        temp=Lbase; Lbase=Rbase; Rbase=temp;
        temp=Ltip; Ltip=Rtip; Rtip=temp;
        temp=Llead; Llead=Rlead; Rlead=temp;
        temp=Ltrail; Ltrail=Rtrail; Rtrail=temp;
    case 3 %digitizing calibration was right-hand
        %all Y values should be the other way around
        head(2,:) = -head(2,:);     %1
        abd(2,:) = -abd(2,:);       %2
        Lbase(2,:) = -Lbase(2,:);   %3
        Ltip(2,:) = -Ltip(2,:);     %4
        Rbase(2,:) = -Rbase(2,:);   %5
        Rtip(2,:) = -Rtip(2,:);     %6
        Llead(2,:) = -Llead(2,:);   %7
        Rlead(2,:) = -Rlead(2,:);   %8
        Ltrail(2,:) = -Ltrail(2,:); %9
        Rtrail(2,:) = -Rtrail(2,:); %10
        %leading edge and trailing edge should be replaced
        temp=Llead; Llead=Ltrail; Ltrail=temp;
        temp=Rlead; Rlead=Rtrail; Rtrail=temp;
end
%% Store original extended coords
coords(1).extended=head'; 
coords(2).extended=abd'; 
coords(3).extended=Lbase'; 
coords(4).extended=Ltip'; 
coords(5).extended=Rbase'; 
coords(6).extended=Rtip'; 
coords(7).extended=Llead'; 
coords(8).extended=Rlead'; 
coords(9).extended=Ltrail'; 
coords(10).extended=Rtrail'; 
%% Scale correction
xyz=((head'+abd')./2); %coords of body center
v=velocity(xyz,Hz);
if mean(sqrt((v(:,1).^2)+(v(:,2).^2)+(v(:,3).^2)))>10 %if mean 3D velocity > 10 m/s, it means that calibration is in cm (in Eyal's data). Few insect's species can achieve this speed (according to Gal)
    for dig_point=1:10
        coords(dig_point).extended=coords(dig_point).extended./100;
    end
    fprintf('Seems that in ''%s'' calibration is in cm instead of meter (3D speed is >10 m/s).\nHence, data is devided by 100.\n',data.save.filename);
end
%% Removing extra data frames 
for dig_point=1:10
    coords(dig_point).original=coords(dig_point).extended((1+prior_frames):(end-posterior_frames),:);
end
%% Filter original data
if use_filter
    coords = filter_digitizing_points(coords,Hz,prior_frames,posterior_frames); %filter % store
    for dig_point=1:size(coords,2) %for each digitizing point
        tmp_coords(dig_point).data=reshape(coords(dig_point).filtered',1,3,size(coords(dig_point).filtered,1)); %reshape to the format in kine's data variable
        coords(dig_point).defualt=coords(dig_point).filtered; %this will used for later data processing
        coords(dig_point).extended((1+prior_frames):(end-posterior_frames),:)=coords(dig_point).filtered; %update extended with filterred coords
    end 
else
    for dig_point=1:size(coords,2) %for each digitizing point
        tmp_coords(dig_point).data=reshape(coords(dig_point).original',1,3,size(coords(dig_point).original,1)); %reshape to the format in kine's data variable
        coords(dig_point).defualt=coords(dig_point).original; %this will used for later data processing
    end 
end %update kine's data variable with the updated coords
data.kine.body.data.coords(1,:,start_frame:end_frame)=tmp_coords(1).data;
data.kine.body.data.coords(2,:,start_frame:end_frame)=tmp_coords(2).data; 
data.kine.left_wing.data.coords(1,:,start_frame:end_frame)=tmp_coords(3).data;
data.kine.left_wing.data.coords(2,:,start_frame:end_frame)=tmp_coords(4).data; 
data.kine.right_wing.data.coords(1,:,start_frame:end_frame)=tmp_coords(5).data;
data.kine.right_wing.data.coords(2,:,start_frame:end_frame)=tmp_coords(6).data; 
data.kine.left_wing.data.model_coords(:,13,start_frame:end_frame)=tmp_coords(7).data;
data.kine.right_wing.data.model_coords(:,13,start_frame:end_frame)=tmp_coords(8).data;
data.kine.left_wing.data.model_coords(:,31,start_frame:end_frame)=tmp_coords(9).data; 
data.kine.right_wing.data.model_coords(:,31,start_frame:end_frame)=tmp_coords(10).data;
%% interpollate cycle's coordinates to uniform data size
for dig_point=1:size(coords,2)
    data_size=size(coords(dig_point).defualt,1);
    t = (1./data_size:1./data_size:1.0)'; %will be used for intepolations: X values (breakpoints) between 0 to 1 with data (Y value) in the wing cycle
    for dim=1:3 %repeat for each axis 
        cs(dig_point).spline(dim)=spline(t',[coords(dig_point).defualt(1,dim) coords(dig_point).defualt(:,dim)' coords(dig_point).defualt(end,dim)]); % piecewise polynomial form of the cubic spline interpolant of cycle_data for later use with the first and last value in cycle_data are used as the endslopes for the cubic spline.
        coords(dig_point).interpolation(dim,:) = ppval(cs(dig_point).spline(dim),tcorr); %Spline Interpolation with Specified Endpoint Slopes
        coords(dig_point).interpolation(dim,1)=coords(dig_point).defualt(1,dim); coords(dig_point).interpolation(dim,end)=coords(dig_point).defualt(end,dim); %ensure that first and last value in cycle_data are used
    end
end
end