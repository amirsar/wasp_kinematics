function [frames, prior_frames, posterior_frames] = kine_cycles_frames_validation(data, frames, wing_cycle, cycles, continuous) 
        error_tester=1; %start testing in order to make sure there are digitizing data in the cycle
        while error_tester
            if ~data.kine.body.data.coords(end,end,frames(wing_cycle))==0 %if first frame isn't empty
                if ~data.kine.body.data.coords(end,end,frames(wing_cycle+1))==0 %if last frame isn't empty
                    error_tester=0; %validation succesful, stop validation
                else
                    frames(wing_cycle+1)=frames(wing_cycle+1)-1; %if last frame doesn't contain data, reduce by one cycle range
                    fprintf('Error!\n Last frame in cycle %d doesn''t exist. Modified from %d to %d.\n',wing_cycle,(frames(wing_cycle+1)+1),frames(wing_cycle+1));
                    if frames(wing_cycle)==frames(wing_cycle+1) %if cycle is empty
                        uiwait(errordlg('Digitizion doesn''t exist for this cycle !'));
                        break
                    end
                end
            else
                frames(wing_cycle)=frames(wing_cycle)+1; %if first frame doesn't contain data, reduce by one cycle range
                fprintf('Error!\n First frame in cycle %d doesn''t exist. Modified from %d to %d.\n',wing_cycle,(frames(wing_cycle)-1),frames(wing_cycle));
                if frames(wing_cycle)==frames(wing_cycle+1) %if cycle is empty
                    uiwait(errordlg('Digitizion doesn''t exist for this cycle !'));
                    break
                end
            end
        end
        %% add extra frames in edges in order to avoid extrapolatin mistakes (velocity & spline functions)
        prior_frames=0;
        posterior_frames=0;
        if (error_tester==0) %if frames are validated
            if wing_cycle==1 || strcmp(continuous,'f') %if it is the first cycle, or if cycles are fragmented
                extra_data=1; %start adding extra frames before first frame
                tested_frame=frames(wing_cycle)-1; %test the previous frame
                if tested_frame==0 %if tested frame is before the frame range of the solution
                    extra_data=0; %stop adding extra frames
                end
                while extra_data && (prior_frames<5) %add up to 5 frames
                    if ~data.kine.left_wing.data.coords(end,end,tested_frame)==0 %if tested frame contain digitizing data on left wing
                        prior_frames=prior_frames+1; %write down how many frames before contain data
                        tested_frame=tested_frame-1; %test the previous frame
                    else %if tested frame don't contain data
                        extra_data=0; %stop adding extra frames
                    end
                end
            else %if it isn't the first cycle
                if ~strcmp(continuous,'f') %if cycles are continuous
                    prior_frames=5; %there should be 5 frames with data in the prior_frames
                else %if cycles are fragmented
                    
                end
            end
            if wing_cycle==cycles || strcmp(continuous,'f') %if it is the last cycle, or if cycles are fragmented
                extra_data=1; %start adding extra frames after last frame
                tested_frame=frames(wing_cycle+1)+1; %test the posterior frame
                if tested_frame>size(data.kine.left_wing.data.coords,3) %if tested frames are beyond the the frame range of the solution
                    extra_data=0; %stop adding extra frames
                end
                while extra_data  && (posterior_frames<5) %add up to 5 frames
                    if ~data.kine.left_wing.data.coords(end,end,tested_frame)==0 %if tested frame contain digitizing data on left wing
                        posterior_frames=posterior_frames+1; %write down how many frames after contain data
                        tested_frame=tested_frame+1; %test the next frame
                    else %if tested frame don't contain data
                        extra_data=0; %stop adding extra frames
                    end
                end
            else %if it isn't the last cycle
                if ~strcmp(continuous,'f') %if cycles are continuous
                    posterior_frames=5; %there should be 5 frames with data in the posteriorframes
                end
            end
        end
end