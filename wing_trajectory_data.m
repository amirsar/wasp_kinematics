%=========================================================================%
% 21/11/2021 By Amir Sarig                                                %
%                                                                         %
% - Calculate & store 3D wing trajectory with body possition and chord.   %
%                                                                         %
% - Make use of functions: wing_trajectory_plot                           %
%                                                                         %
% Called from 'AnglesExcel_Amir' & 'extract_data_wasp_forward_flight'     %
% Calls 'wing_trajectory_plot' func.                                      %
%=========================================================================%

function [rotated_coords,coords_index,stroke_plane]=wing_trajectory_data(coords,ind,idx,rotated_coords,coords_index,draw_plot,tcorr,V_tip_storage,stroke_plane,V_data,elevation)
%% Obtaining body & wing points
head = coords(1).defualt';
abd = coords(2).defualt';
Lbase = coords(3).defualt';
Ltip = coords(4).defualt';
Rbase = coords(5).defualt';
Rtip = coords(6).defualt';
Llead = coords(7).defualt';
Rlead = coords(8).defualt';
Ltrail = coords(9).defualt';
Rtrail = coords(10).defualt';
%% Calculation
% Defining the axes of the body
x = (head-abd)./repmat(sqrt((head(1,:)-abd(1,:)).^2+(head(2,:)-abd(2,:)).^2+(head(3,:)-abd(3,:)).^2),3,1); %body length unit vector
yy = (Lbase-Rbase)./repmat(sqrt((Lbase(1,:)-Rbase(1,:)).^2+(Lbase(2,:)-Rbase(2,:)).^2+(Lbase(3,:)-Rbase(3,:)).^2),3,1); %body width unit vector
tet = acos(dot(yy,x)); %angle between x & yy (both are already unit vectors)
zz = cross(x,yy); %body hight vector
z = zz./repmat(sin(tet),3,1);  %body hight unit vector (xXyy=1*1*sin(tet))
y = cross(z,x);  %corrected body width unit vector (so it is orthogonal to x)
TM = permute(reshape([x;y;z],3,3,length(head)),[2,1,3]);  %rotation matrix (each individual 2D matrix is [x1,x2,x3;y1,y2,y3;z1,z2,z3])
AllD = reshape([head;abd;Lbase;Ltip;Rbase;Rtip;Llead;Rlead;Ltrail;Rtrail],3,10,length(head));  %All Data
AllD = AllD-repmat((AllD(:,1,:)+AllD(:,2,:))/2,1,10); %normalize by reducing body center from all digitizing points
% Rotation
P = zeros(size(AllD));
for n = 1:length(head)
    P(:,:,n) = TM(:,:,n)*AllD(:,:,n);
end
%data interpolation
if coords_index(1,1)==0 %how many rows are in 'coords_index' +1
    coords_row=1;
else
    coords_row=size(coords_index,1)+1; %wing cycle number
end
for dig_pt=1:size(AllD,2) %repeat interpolation for each calculated point
    cycle_data=reshape(P(:,dig_pt,:),3,length(head)); %extract left wing tip data
    [r,rr] = size(cycle_data);
    t = (1./rr:1./rr:1.0)'; %X values (breakpoints) between 0 to 1 with data (Y value) in the wing cycle
    for i=1:r %repeat for each axis (r should always be equal 3)
        cs(i)=spline(t,[cycle_data(i,1) cycle_data(i,:) cycle_data(i,end)]); % piecewise polynomial form of the cubic spline interpolant of cycle_data for later use with the first and last value in cycle_data are used as the endslopes for the cubic spline.
        corrected_data(i,:) = ppval(cs(i),tcorr); %Spline Interpolation with Specified Endpoint Slopes
        corrected_data(i,1)=cycle_data(i,1); corrected_data(i,end)=cycle_data(i,end); %ensure that first and last value in cycle_data are used
    end
    %store data one cycle after the other. order: head;abd;Lbase;Ltip;Rbase;Rtip;Llead;Rlead;Ltrail;Rtrail
    rotated_coords.x_axis(dig_pt).data(coords_row,:)=corrected_data(1,:); %store X coordinates along cycle
    rotated_coords.y_axis(dig_pt).data(coords_row,:)=corrected_data(2,:); %store Y coordinates along cycle
    rotated_coords.z_axis(dig_pt).data(coords_row,:)=corrected_data(3,:); %store Z coordinates along cycle
end %end of interpolation of one digitizing point
coords_index(coords_row,:)=[ind,idx]; %store the indentity of each coordinates line

%% stroke plain
body_velocity = [mean(x,2)';mean(y,2)';mean(z,2)']*[V_data{coords_row,7};V_data{coords_row,8};V_data{coords_row,9}]; %rotate 3D velocity into body oriantation
for i=2:3 %for each wing (left and right)
    wing_tip=reshape(P(:,i*2,:),3,length(head)); %extract cycle's wing tip data
    [maxZ(i-1),maxZ_index(i-1)]=max(wing_tip(3,:)); %max Z value and it's index in this cycle
    [minZ(i-1),minZ_index(i-1)]=min(wing_tip(3,:)); %min Z value and it's index in this cycle
    stroke_XZ_plane((i-1),:)=[wing_tip(1,maxZ_index(i-1))-wing_tip(1,minZ_index(i-1)),maxZ(i-1)-minZ(i-1)]; %max(x,z)-min(x,z)
    %wing_length_vectors = reshape(P(:,i*2,:)-P(:,(i*2)-1,:),3,length(kine_frames));  %wing length vectors (tip-base) of all the wing cycle
    %amplitude(i-1) = atan2d(wing_length_vectors(3,maxZ_index(i-1)),-wing_length_vectors(2,maxZ_index(i-1)))-atan2d(wing_length_vectors(3,minZ_index(i-1)),-wing_length_vectors(2,minZ_index(i-1))); %amplitude
    amplitude(i-1)=radtodeg(max(elevation(i-1,:))-min(elevation(i-1,:))); %amplitude
end
mean_stroke_plane=mean(stroke_XZ_plane,1); %right and left wings mean stroke plane vector XZ
stroke_plane_angle_from_velocity=radtodeg(acos(dot(body_velocity([1 3])',mean_stroke_plane)./(norm(body_velocity([1 3])')*norm(mean_stroke_plane)))); %stroke plane angle from velocity vector: teta=acos((A dot B)/(norm(A)*norm(B)))
stroke_plane_angle_from_body=radtodeg(acos(dot([mean(x(1,:)) mean(x(3,:))],mean_stroke_plane)/(norm([mean(x(1,:)) mean(x(3,:))])*norm(mean_stroke_plane)))); %stroke plane angle from body vector: teta=acos((A dot B)/(norm(A)*norm(B)))
stroke_plane{coords_row,:}=[V_data{coords_row,1},V_data{coords_row,2},wing_tip(1,maxZ_index(1)),maxZ(1),wing_tip(1,minZ_index(1)),minZ(1),wing_tip(1,maxZ_index(2)),maxZ(2),wing_tip(1,minZ_index(2)),minZ(2),amplitude(1),amplitude(2),mean_stroke_plane(1),mean_stroke_plane(2),stroke_plane_angle_from_velocity,stroke_plane_angle_from_body,body_velocity(1),body_velocity(2),body_velocity(3)]; %store data
if draw_plot %if it is the last cycle
    stroke_plane.Properties.VariableNames={'Trial','Cycle','LmaxX','LmaxZ','LminX','LminZ','RmaxX','RmaxZ','RminX','RminZ','Lamplitude','Ramplitude','Avg_SP_vector_X','Avg_SP_vector_Z','SP_angle2velocity','SP_angle2body','rotatedVbody_X','rotatedVbody_Y','rotatedVbody_Z'}; %add labels to data
end
%% draw wingtip trajectory
if draw_plot %if enter here, it means that all data for figure exist in 'coords'
    trials=unique(coords_index(:,1)); %exctract amount of trials
    for index=1:length(trials) %for each trial
        for dig_pt=1:length(rotated_coords.x_axis)  %calc & store each trial mean of each digitizing point in each time-point. mean is calculated along columns
            mean_x(dig_pt).data(index,:)=mean(rotated_coords.x_axis(dig_pt).data(coords_index(:,1)==trials(index),:),1);
            mean_y(dig_pt).data(index,:)=mean(rotated_coords.y_axis(dig_pt).data(coords_index(:,1)==trials(index),:),1);
            mean_z(dig_pt).data(index,:)=mean(rotated_coords.z_axis(dig_pt).data(coords_index(:,1)==trials(index),:),1);
        end
    end
    wing_trajectory_plot(tcorr,mean_x,mean_y,mean_z,V_tip_storage) %plot mean wingtip trajectory
    t=sprintf('Wing tip trajectory based on %d flapping cycles from %d individuals',length(coords_index),index);
    title(t) %add title to figure
end
end