%% load excel data so can use to locate tif files
%% add filepaths, initialize data that we are looking for

clear

load('excelshapedata.mat')

%% initialize ratmap head

load RatMap_Head;
headx = RatMap_Head.x'; %x-coordinates of original rat head, row vector due to transpose
heady = RatMap_Head.y';
headz = RatMap_Head.z';
head_data = [headx;heady;headz];

%bring ratmaphead data to have nose (or a point lil off the nose since we
%are corresponding to tracked nose pt which is never ON the nose per se) on
%the origin. = translate nose to origin
t_mat=eye(4);
t_mat(1:3,4) = [0;-1.3993;0]; %nose is at [0,1.3993,0] -- find using max(heady) -- so want to translate so this pt, and all attached pts move to origin
head_data_pre_t = [head_data;ones(1,size(head_data,2))];
t_head_data=t_mat*head_data_pre_t;
t_head_data(4,:) = [];

%pitch head to have pitch equal to zero
%nose is now at 0,0,0

%better pitch reference (bridge of head/nose NOT the eye mdpt.
bridgeofhead_pt = [0,-20,2.8];
t_ratmap_pitch_ref = t_mat*[bridgeofhead_pt';1];
t_ratmap_pitch_ref(4,:) = [];

%zero_pitch = -1*atan2d(t_ratmap_eye_mdpt(3),t_ratmap_eye_mdpt(2))- 180;
zero_pitch = 1*atan2d(t_ratmap_pitch_ref(3),-t_ratmap_pitch_ref(2));
zeropitch_r_mat = getRX(2*pi*zero_pitch/360);
zeropitch_ratmap_head = zeropitch_r_mat*t_head_data;

%rotate head to default 0 degree position
rot_90deg= getRZ(2*pi*90/360);
initialized_head_data = rot_90deg*zeropitch_ratmap_head;

%% make and save videos with tif + tracked pts + head angle plotted line

load('td_headpoints') %top down view
load('tracked_headpitch') %side view

close all

figure(1);
f1 = figure(1);
set(f1, 'Units', 'normalized', 'Position', [0,0,1,1]);

for clip_counter=1:length(td_headpoints)

    clip_data = td_headpoints(clip_counter).trackedPts;
    
    
    %%create pitch array or this clip
    
    %%for right now just dealing with clip 1, eventually want to write
    %%loops to make sure you only plot the clips for which you have both
    %%headpoints in top view(yaw) and headpoints in side view
    %%(pitch)
    
    %%complete interpolated pitch array for each frame will be contained in
    %%"interpolated_pitch_for_clip"
    interpolated_pitch_for_clip = [];
    tracked_pitch_for_clip = headpitch(clip_counter).vheadang;
    pitch_frame_num = headpitch(clip_counter).frameNum;
    for n=1:length(pitch_frame_num);
        %assumes starting pitch is at frame #1, otherwise will have to
        %rewrite interpolation, still easy)
        if n==length(pitch_frame_num)
            n_frames_in_clip = size(clip_data,3);
            piecewise_pitch = linspace(tracked_pitch_for_clip(n), 0,n_frames_in_clip - pitch_frame_num(n)+1);
            piecewise_pitch = piecewise_pitch(2:end);
        else
            piecewise_pitch = linspace(tracked_pitch_for_clip(n),tracked_pitch_for_clip(n+1),pitch_frame_num(n+1) - pitch_frame_num(n)+1);
            piecewise_pitch = piecewise_pitch(2:end-1);
        end
        interpolated_pitch_for_clip= [interpolated_pitch_for_clip;tracked_pitch_for_clip(n);piecewise_pitch'];
    end
    
    %%interpolate tracked points for side view
    
    %%(only need to interpolate one point of the two (nose and eye) because
    %%the interpolation of tracked pts may not correspond to the
    %%interpolation of head pitches, so we are choosing the NOSE point to
    %%be interpolated, since it makes it easier to correspoond to our
    %%ratmaphead which has its nose point tracked during various
    %%rotations/scaling/translations, etc.
    
    side_clip_data = headpitch(clip_counter).trackedPts;%%%again will need to make sure you are accessing the correct indices for side and topdown data respectively
    interpolated_side_clip_data = [];
        for n=1:length(pitch_frame_num)
        %assumes starting pitch is at frame #1, otherwise will have to
        %rewrite interpolation, still easy).  ADDITIONALLY, this assumes
        %last index is last tracked pt, rather than for head pitch where we
        %interpolate for ALL frames, not just the # frames that we have
        %side view-tracked pts for.
            if n==length(pitch_frame_num)
                interpolated_side_clip_data= [interpolated_side_clip_data;side_clip_data(1,:,n)];
            else
                piecewise_trackednoseX_side_view = linspace(side_clip_data(1,1,n),side_clip_data(1,1,n+1),pitch_frame_num(n+1) - pitch_frame_num(n)+1);
                piecewise_trackednoseY_side_view = linspace(side_clip_data(1,2,n),side_clip_data(1,2,n+1),pitch_frame_num(n+1) - pitch_frame_num(n)+1);
                interpolated_side_clip_data= [interpolated_side_clip_data;side_clip_data(1,:,n);[piecewise_trackednoseX_side_view(2:end-1); piecewise_trackednoseY_side_view(2:end-1)]'];
            end
%             pause;
        end

        %find average distance between eye and nose in the tracked frames
        %for side view.  This will be used to plot a "pitch line" in the
        %side views
        dist = [];
        for n = 1:size(side_clip_data,3)
            dist = [dist, sqrt(sum((side_clip_data(1,:,n) - side_clip_data(2,:,n)).^2))];
        end
        avg_side_dist = sum(dist)/size(side_clip_data,3);
        
        %%%%%%%%%%Calculate scaling factor for top down view for THIS CLIP
        
        %%take scaling from frame with least head pitch
        min_pitch = min(abs(interpolated_pitch_for_clip));
        i_min_pitch= find(abs(interpolated_pitch_for_clip)==min_pitch);
        disp(strcat('The frame with the least pitch is frame # ', num2str(i_min_pitch), ' with pitch equal to ', num2str(interpolated_pitch_for_clip(i_min_pitch))))
        midear = [(clip_data(1,1,i_min_pitch) + clip_data(3,1,i_min_pitch))/2, (clip_data(1,2,i_min_pitch) + clip_data(3,2,i_min_pitch))/2];
        nose = clip_data(2,:,i_min_pitch);
        td_scale_for_clip = sqrt(sum((midear-nose).^2));
        clear i_min_pitch
        
    for fcounter = 1:size(clip_data,3)
        
        %change working directory to correct tif folder corresponding to trial
        %number
        
        trial_num = td_headpoints(clip_counter).trialNum;
        tif_cell = strcat({'rat '}, {num2str(Rat(trial_num))}, {' trial '}, num2str(Trl(trial_num)), {' '},char(Perf(trial_num+2)));
        tiffolder = char(tif_cell);
        
        if length(num2str(Date(trial_num)))==5
            trial_date = strcat('0',num2str(Date(trial_num)));
        elseif length(num2str(Date(trial_num)))==6
            trial_date = num2str(Date(trial_num));
        end
        
        cd(strcat('F:\ShapeData_2015_Working\TIF_files\',trial_date,'\',tiffolder))
        
        %plot tif
        
        subplot(2,3,1)
        
        tif_frame = td_headpoints(clip_counter).startFrame+fcounter-1;
        tif_prefix = 'c001s0001';
        tifname = [tif_prefix int2str2(tif_frame,100000)];
        img = loadtif(tifname);
        imagesc(flipud(fliplr((img)))); colormap('gray');
        set(gca,'YDir','reverse')
        axis square
        ho;
        
        %calculate angle of head based on tracked points
        
        ax= gca;
        
        if ax.YDir == 'reverse'
            nose_pt = clip_data(2,:,fcounter);
            eye_mdpt = [(clip_data(1,1,fcounter) + clip_data(3,1,fcounter))/2, (clip_data(1,2,fcounter) + clip_data(3,2,fcounter))/2];
            eye_mdpt_t = eye_mdpt-nose_pt;
            rot_rad = atan2d(eye_mdpt_t(2),-eye_mdpt_t(1));
        elseif ax.YDir=='normal'
            nose_pt = clip_data(2,:,fcounter);
            eye_mdpt = [(clip_data(1,1,fcounter) + clip_data(3,1,fcounter))/2, (clip_data(1,2,fcounter) + clip_data(3,2,fcounter))/2];
            eye_mdpt_t = eye_mdpt-nose_pt;
            rot_rad = atan2d(eye_mdpt_t(2),eye_mdpt_t(1)); %atan2(y' coord, x'coord') = 4 quadrant inverse tangent
        end
        

        %plot tracked points
        
        plot(clip_data(1,1,fcounter),clip_data(1,2,fcounter), 'c*');ln1
        plot(clip_data(2,1,fcounter),clip_data(2,2,fcounter), 'y*');ln1
        plot(clip_data(3,1,fcounter),clip_data(3,2,fcounter), 'g*');ln1
        ho;
        
        %plot head angle line
        plot([nose_pt(1),eye_mdpt(1)],[nose_pt(2),eye_mdpt(2)],'y-');ln2
        text(150,50,char(strcat({'Rotation = '},num2str(rot_rad, '%.2f'), {' degrees'})),'color','white')
        title(char(strcat({'Trial No. '}, num2str(td_headpoints(clip_counter).trialNum), {' Frame No. '}, num2str(fcounter))), 'FontSize',10)
        xlabel('x');ylabel('y');
        hold off;
        
        xlim_for_fig2 = ax.XLim;
        ylim_for_fig2 = ax.YLim;

        %%%%%%%%%%%%%%%%%%%%%% plot ratmap head in correct positions + tracked pts
        
        %scale ratmaphead
        
        %%%these below position pts correspond to the
        %%%'zeropitch_ratmap_head' (unrotated).  Therefore, they do not
        %%%corresond to the initialized ratmap at the start of this program
        %%%titled "initialized_head_data"
        ratmap_rightear_pt = [19.1,-49,-6.5];
        ratmap_leftear_pt= [-16.4,-47.5,-6.5];
        %plot these over zeropitch_ratmap_head in order to see what points
        %I am using to the ultimately correspond to the tracked points.
        %(this affects scaling)
        
        %the z coordinates of these ratmap points^^ is not exactly relevant, because we are
        %only using the X/Y coordinates to scale for the top down view anyways
        ratmap_ear_mdpt = (ratmap_leftear_pt + ratmap_rightear_pt)/2;
        ratmap_nose_pt = [0,.4228,0]; %find using max(zeropitch_ratmap_head(2,:))
        ratmap_midline_distance = sqrt((ratmap_nose_pt(1) - ratmap_ear_mdpt(1))^2+(ratmap_nose_pt(2)-ratmap_ear_mdpt(2))^2);
        
        %this scale factor is constant for each frame
        scale_factor = td_scale_for_clip/ratmap_midline_distance;
        scaled_head = scale_factor*initialized_head_data;
        
        %pitch scaled head
        %%%this time the pitch head will rotate about the y axis and not the x-axis, since we
        %%%have positione the ratmap in a different resting position that
        %%%what teh default ratmaphead is positioned as.
        
        %%pitch the head by the experimental value
        exp_pitch_deg = 1*interpolated_pitch_for_clip(fcounter);
        %-1 multiplied above^^^^ because these pitch angles correspond to the pitch
        %between the eye and the nose, and we are rotating about the nose,
        %so we actually want the opposite direction of rotation in order for the correct pitch to be applied
        r_mat_exp_pitch = getRY(2*pi*exp_pitch_deg/360);
        exp_pitch_ratmap_head = r_mat_exp_pitch*scaled_head;    
        
        %rotate scaled, pitched ratmaphead to rotation angle for this frame
        rot_rad2 = atan2(eye_mdpt_t(2),eye_mdpt_t(1));
        r_mat_atframe = getRZ(rot_rad2);
        r_scaled_head_atframe = r_mat_atframe*exp_pitch_ratmap_head;
        
        %translate scaled head to have the scaled nose with the tracked
        %nose pt
        r_scaled_head_atframe_pre_t = [r_scaled_head_atframe;ones(1,size(r_scaled_head_atframe,2))];
        t_mat_atframe = eye(4);
        t_mat_atframe(1:3,4) = [nose_pt,0]';
        t_r_scaled_head_atframe = t_mat_atframe*r_scaled_head_atframe_pre_t;
        
        %plot pitched, scaled, rotated, and translated ratmap head

        subplot(2,3,2)
        
        h = trisurf(RatMap_Head.t,t_r_scaled_head_atframe(1,:),t_r_scaled_head_atframe(2,:),t_r_scaled_head_atframe(3,:), 'FaceColor','k','EdgeColor','none','FaceAlpha',.1);
        ho;
        view(2)
        xlim(xlim_for_fig2);
        ylim(ylim_for_fig2);
        ax2 = gca;
        ax2.YDir = 'reverse';
        axis square
        
        %plot tracked points
        plot(clip_data(1,1,fcounter),clip_data(1,2,fcounter), 'c*');ln1
        plot(clip_data(2,1,fcounter),clip_data(2,2,fcounter), 'y*');ln1
        plot(clip_data(3,1,fcounter),clip_data(3,2,fcounter), 'g*');ln1
        title(char(strcat({'Trial No. '}, num2str(td_headpoints(clip_counter).trialNum), {' Frame No. '}, num2str(fcounter), {' with ratmap head'})), 'FontSize',10)
        xlabel('x');ylabel('y');
        hold off;

%%%%%%%%%%%%%%%%%%%plot overlay fo ratmap and tif file

        % plot tif image
        
        subplot(2,3,3)
        img = loadtif(tifname);
        imagesc(flipud(fliplr((img)))); colormap('gray');
        ho;


        %plot tracked points
        ho;
        plot(clip_data(1,1,fcounter),clip_data(1,2,fcounter), 'c*');ln1
        plot(clip_data(2,1,fcounter),clip_data(2,2,fcounter), 'y*');ln1
        plot(clip_data(3,1,fcounter),clip_data(3,2,fcounter), 'g*');ln1
        ho;
        
        % plot ratmap overlay
        %%%%%plotting downsampled ratmap head to allow to see through to
        %%%%%tif image
        h = plot(t_r_scaled_head_atframe(1,1:5:end),t_r_scaled_head_atframe(2,1:5:end),'.', 'MarkerSize',3);
        ho;
        xlim(xlim_for_fig2);
        ylim(ylim_for_fig2);
        ax2 = gca;
        ax2.YDir = 'reverse';
        axis square
        
        title(char(strcat({'Trial No. '}, num2str(td_headpoints(clip_counter).trialNum), {' Frame No. '}, num2str(fcounter), {' overlay of ratmap head'})), 'FontSize',10)
        xlabel('x');ylabel('y');
        hold off;

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot ratmaphead in side view



%%%%%%%%%%%%plot side view sideview tif with tracked pts

        subplot(2,3,4)
        
        %plot tif
        tif_prefix = 'c002s0001'; %camera c002, side view
        tifname_side = [tif_prefix int2str2(tif_frame,100000)];
        img_side = loadtif(tifname_side);
        imagesc(img_side); colormap('gray');
        %%%%%%%%%%%NOTICE THE DIFFERENT AXIS LABELING, IT'S MEANT TO MATCH
        %%%%%%%%%%%the ratmap axes which actually are showing the y and z
        %%%%%%%%%%%dimensions, so taht you can compare from the top down
        %%%%%%%%%%%(x/y) to side view (y/z) dimensions
        xlabel('y');ylabel('z');
        ho;

        %plot tracked points
        plot(interpolated_side_clip_data(fcounter,1),interpolated_side_clip_data(fcounter,2), 'y*');ln1%nose
        ho;
        
        %plot pitch line
        zero_pitch_avg_line = [linspace(0,avg_side_dist,100);zeros(1,100);zeros(1,100)];%third zero vec is just to make 3D so rotation multiplication works
        rot_pitch_line_mat = getRZ(2*pi*-1*interpolated_pitch_for_clip(fcounter)/360);
        rot_pitch_line = rot_pitch_line_mat*zero_pitch_avg_line;       
        t_mat = eye(4);
        t_mat(1:2,4) = interpolated_side_clip_data(fcounter,:)' - rot_pitch_line(1:2,end);
        rot_pitch_line_pre_t = [rot_pitch_line;ones(1,size(rot_pitch_line,2))];
        final_pitch_line = t_mat*rot_pitch_line_pre_t;
        final_pitch_line(4,:) = [];
        
        plot(final_pitch_line(1,:),final_pitch_line(2,:),'y');ln4;
        
        ax_row2 = gca;

%%%%%%%%%%%%plot ratmaphead in side view
        
        %scale ratmaphead
        
%         %%%these below position pts correspond to the default ratmap when simply plotting
%         %%%ratmap head using Plot_RatMap_Head.  They do not corresond to
%         %%%the initialized ratmap at the start of this program titled
%         %%%"initialized_head_data"
%         ratmap_lefteye_pt = [-14.9,-32,-10];
%         ratmap_righteye_pt= [15.5,-32,-10];
%         ratmap_eye_mdpt = (ratmap_lefteye_pt + ratmap_righteye_pt)/2;
%         ratmap_nose_pt = [0,2,0];
%         %%%NOTE: THE BELOW DISTANCES ARE DIFFERENT THAN THE SCALING USED
%         %%%FOR THE TOP DOWN VIEW, specifically, the Y and Z coordinates are
%         %%%used for side view, whereas X and Y distances are used for top
%         %%%down.  this is because we have orthoganal views of cameras,
%         %%%therefoe orthogonal (and different therefore) dimensions that we
%         %%%are scaling accoding to.
%         ratmap_midline_distance = sqrt((ratmap_nose_pt(3) - ratmap_eye_mdpt(3))^2+(ratmap_nose_pt(2)-ratmap_eye_mdpt(2))^2);
%         tracked_midline_distance = avg_side_dist;
%         scale_factor = tracked_midline_distance/ratmap_midline_distance;
%         scaled_head = scale_factor*initialized_head_data;

        %easier scalng head (by just glancing, aproximate)
        scale_factor = 4.4;
        scaled_head = scale_factor*initialized_head_data;
        

        %pitch scaled head
        %%%this time the pitch head will rotate about the y axis and not the x-axis, since we
        %%%have positione the ratmap in a different resting position that
        %%%what teh default ratmaphead is positioned as.
        
        %%pitch the head by the experimental value
        exp_pitch_deg = 1*interpolated_pitch_for_clip(fcounter);        
        %-1 multiplied above^^^^ because these pitch angles correspond to the pitch
        %between the eye and the nose, and we are rotating about the nose,
        %so we actually want the opposite direction of rotation in order for the correct pitch to be applied
        r_mat_exp_pitch = getRY(2*pi*exp_pitch_deg/360);
        exp_pitch_ratmap_head = r_mat_exp_pitch*scaled_head;    

        %rotate 180 degrees about y axis
        roty = getRY(2*pi*180/360);
        exp_pitch_ratmap_head = roty*exp_pitch_ratmap_head;
%         
%         %rotate 180 degrees about x axis
%         rotx = getRX(2*pi*180/360);
%         exp_pitch_ratmap_head = rotx*exp_pitch_ratmap_head;

        %rotate scaled, pitched ratmaphead to rotation angle for this frame
        rot_rad2 = atan2(eye_mdpt_t(2),eye_mdpt_t(1));
        r_mat_atframe = getRZ(rot_rad2);
        r_scaled_head_atframe = r_mat_atframe*exp_pitch_ratmap_head;
        
        %translate scaled head to have the scaled nose with the tracked nose pt
        r_scaled_head_atframe_pre_t = [r_scaled_head_atframe;ones(1,size(r_scaled_head_atframe,2))];
        t_mat_atframe = eye(4);
        t_mat_atframe(1:3,4) = [0,interpolated_side_clip_data(fcounter,1),interpolated_side_clip_data(fcounter,2)]';
        t_r_scaled_head_atframe = t_mat_atframe*r_scaled_head_atframe_pre_t;
        
        subplot(2,3,5)
        
        %plot ratmaphead
        h_side = trisurf(RatMap_Head.t,t_r_scaled_head_atframe(1,:),t_r_scaled_head_atframe(2,:),t_r_scaled_head_atframe(3,:), 'FaceColor','k','EdgeColor','none','FaceAlpha',.1);
        axis square
        ylim(ax_row2.XLim);zlim(ax_row2.YLim);

        title(char(strcat({'Trial No. '}, num2str(td_headpoints(clip_counter).trialNum), {' Frame No. '}, num2str(fcounter), {' with ratmap head'})), 'FontSize',10)
        xlabel('x');ylabel('y');zlabel('z')
        a2 = gca;
        a2.XDir = 'reverse';
        a2.ZDir = 'reverse';
        view(90,0)
        
        
        %plot tracked pt
        ho;
        plot3(0,interpolated_side_clip_data(fcounter,1),interpolated_side_clip_data(fcounter,2), 'y*');ln1%nose
        hold off;


%%%%%%%%%%%%%overlay side view tif and side view ratmap

    subplot(2,3,6)
    
    
        %plot tif
        img_side = loadtif(tifname_side);
        imagesc(img_side); colormap('gray');
        ylabel('z');xlabel('y');
        ho;
        
        %plot ratmaphead
        h_side = plot(t_r_scaled_head_atframe(2,1:5:end),t_r_scaled_head_atframe(3,1:5:end),'.', 'MarkerSize',3);
        axis square
        title(char(strcat({'Trial No. '}, num2str(td_headpoints(clip_counter).trialNum), {' Frame No. '}, num2str(fcounter), {' with ratmap head'})), 'FontSize',10)
        
        %plot tracked pt
        ho;
        plot(interpolated_side_clip_data(fcounter,1),interpolated_side_clip_data(fcounter,2), 'y*');ln1%nose
        hold off;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%capture frame
        M3(fcounter) = getframe(gcf);

        if fcounter == length(interpolated_side_clip_data)
            disp(num2str(fcounter))
            break;
        end
    end
    
    %save video
    cd('C:\Users\Sagar P\Documents\hartmann\ratmap_input_from_behav_data\top_down_ratmap_head_comparison_and_original_videos')
    c=char(strcat({'Trial Number '}, num2str(td_headpoints(clip_counter).trialNum),{' both camera views corr'}));
    v = VideoWriter(c);
    v.FrameRate = 5;
    open(v)
    writeVideo(v,M3)
    close(v)
    
    clear M3
   
    if clip_counter==1
        break;
    end
    
end



