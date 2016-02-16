% RatHeadTrack.m
%     A script to track the head and whiskers of a rat from a series of
%     .tiff files. This was coded for Winter2016 BME462.
%   
%   Since this is a collaborative project, please remember to keep
%   everything well commented (womp womp). Current todos will be tracked at
%   the beginning of the file, and if you make any major changes add your
%   initials next to it.


%   Current ToDos:
%       * Generalize naming in photo read-in
%       * Possibly improve edge detection?
%       * Identify whiskers and head
%       * Loading in RatMap based rat locations
%       * 
function RatHeadTrackingSideView
    clc
    close all;
    %Genaralized photo read from directory (create a folder 'video' that stores all
    %videos)
    path = [pwd '/video/*.tif'];
    files = dir(path);
    neighborhood=10;

    % Prepare the new file.
    vidObj = VideoWriter('test.avi');
    open(vidObj);

    % Create an animation.
    set(gca,'nextplot','replacechildren');
    
    global nose_coord head_coord eye_coord
    nose_coord = [];
    head_coord = [];
    eye_coord = [];
    all_radii=[];
    
    tol = 25;
    
    for i = 1:length(files)
    % Runs through all of the pictures for a specific view, converts them to
    % black and white, edge detection, then converts them all into a video.     

        %if the string name contains the second video set break from for
        %loop -- UM
        if ~isempty(strfind(files(i).name,'_c002')) 
            break;
        end

        % Reads in the file (arrangement assumption is that all video are
        % stored another a video folder residing in the path that the code
        % is in)
        image_path_1 = [pwd '/video/' files(i).name];
        I=imread(image_path_1);
        image_path_2 = [pwd '/video/' files(i+1).name];
        I_next = imread(image_path_2);
    
        % store current and next image --UM
        [ED,BW] = getEdges(I);        
        [D_ED,D_BW]=getEdges(I_next);
        BW=~BW; D_BW=~D_BW; %invert image
        tracking=figure(1); imshow(BW);
        s=regionprops('table',BW,'centroid',...
            'MajorAxisLength','MinorAxisLength');
        centers=s.Centroid;
        diameters=mean([s.MajorAxisLength s.MinorAxisLength],2);
        radii=diameters/2; 


        %identify nose and head points --UM
         if i==1
             fig = figure(2); h = imshow(ED); hold on 
             [x_start,y_start] = getpts(fig); % hand pick nose and head coordinates
             %(double click ends selection process)
             plot(x_start,y_start,'g'), hold off;
             [nrows,ncols]=size(get(h,'cdata'));
             xdata=get(h,'xdata'); ydata=get(h,'ydata');
             %Convert axes to pixel coordinates for accurate mapping. Was
             %not essential to the given data set
             x_nose = axes2pix(ncols,xdata,x_start(1)); y_nose = axes2pix(nrows,ydata,y_start(1));
             x_head = axes2pix(ncols,xdata,x_start(2)); y_head=axes2pix(nrows,ydata,y_start(2));
             nose=[x_nose y_nose]; head=[x_head,y_head]; nose_head=[nose;head];
             nose_coord = [nose_coord;nose]; head_coord=[head_coord;head];
             L = pdist(nose_head,'euclidean');
         end
         
        %find the edges in the next image
        [r,c] = find(D_ED==1); %find pixel locations of the edge
        E = [c r]; %imread is the opposite coordinates, flipped col(x) and row(y)
        
        %find the distance between the nose point to the edges of the next
        %image.
        for k=1:size(E,1)
            N=[nose;E(k,:)]; H=[head;E(k,:)]; 
            d_nose(k)=pdist(N,'euclidean');
            d_head(k)=pdist(H,'euclidean');
        end
        
        %Find the next nose and head point by finding the closest edge point to the
        %previous coordinates. The first is user specified, the rest is
        %automatic. --UM
        %Issues:
            % * weirdness at the end
            % * not robust to noise
        %Solutions:
            % * a better metric than min dist could be used
            % * possibly using the initial distance between the head point and the nose point
       %Fixed:
            % * better nose and head tracking looking at neighboring pixels
            
      % * improved neighboring pixels based on connectivity in the region
      %improved nose and head tracking
        
        [~,ind_nose]=min(d_nose); [~,ind_head]=min(d_head);
        
        nose_region=D_ED(E(ind_nose,2)-neighborhood:E(ind_nose,2)+neighborhood,...
            E(ind_nose,1)-neighborhood:E(ind_nose,1)+neighborhood);
        head_region=D_ED(E(ind_head,2)-neighborhood:E(ind_head,2)+neighborhood,...
            E(ind_head,1)-neighborhood:E(ind_head,1)+neighborhood);
        [r_nr,c_nr]=find(nose_region==1); [r_hr,c_hr]=find(head_region==1);
        E_NR=[c_nr,r_nr]-neighborhood-1; E_HR=[c_hr,r_hr]-neighborhood-1;
        idx_n=[E_NR(:,1)+E(ind_nose,1), E_NR(:,2)+E(ind_nose,2)];
        idx_h=[E_HR(:,1)+E(ind_head,1), E_HR(:,2)+E(ind_head,2)];

        for j=1:size(idx_h,1)
            temp_h=[idx_h(j,1) idx_h(j,2)];
            for l=1:size(idx_n,1)
                temp_n=[idx_n(l,1) idx_n(l,2)];
                NH = [temp_h;temp_n];
                d(l,j)= pdist(NH,'euclidean');
            end
        end
      
        dum=abs(d-L); [M,rows]=min(dum);[~,min_col]=min(M); min_row=rows(min_col);
        nose=[idx_n(min_row,1) idx_n(min_row,2)]; nose_coord =[nose_coord;nose];
        head=[idx_h(min_col,1) idx_h(min_col,2)]; head_coord=[head_coord;head];
        nose_head=[nose;head]; L=pdist(nose_head,'euclidean');
    
        % find eye using regionprops:
            % * works well in the beginning when the head turns it loses it
            %   a solution I used was to calculate the change in the head
            %   coordinate and use that difference to estimate where it
            %   would be. works terribly, a better way would be to use a
            %   Kalman filter and use the previous run throughs that work
            %   to estimate current eye location
            
        if isempty(eye_coord) 
            radii_temp=sort(radii,'descend');
            eye_radii=radii(radii==radii_temp(2));
            eye_center=centers(radii==radii_temp(2),:);
            figure(tracking);
            viscircles(eye_center,eye_radii);
        else
            dist=[];
            for c=1:length(centers)
                A=[eye_coord(1,:);centers(c,:)];
                dist(c)=pdist(A,'euclidean');
            end
            [val,ind_eye]=min(dist);
            if val<tol %if the center value is within some tolerance distance value, assume that is the eye center
            eye_center=centers(ind_eye,:);
            eye_radii=radii(ind_eye);

            else
            warning('Eye center with region props lost, offset from the head are going to be used for estimation');
            % Could be changed to Kalman Filter for better estimation
            head_dif=head_coord(i-1,:)-head_coord(i,:);
            eye_center=[eye_coord(1,1)-head_dif(1), eye_coord(1,2)-head_dif(2)];
            eye_radii=max(all_radii);
            end
        end
               
        eye_coord=[eye_center;eye_coord];
        all_radii=[eye_radii;all_radii];
        
        %Display nose and head points on image
        figure(tracking); hold on;
        viscircles(eye_center,eye_radii);
        plot(nose_coord(i,1), nose_coord(i,2),'r*',... 
            head_coord(i,1),head_coord(i,2),'go');
        hold off

        frame = getframe;
        writeVideo(vidObj,frame);
        
        dbstop if error %debug mode in case of error
        clearvars d l j dum  
    end

    frame = getframe;
    writeVideo(vidObj,frame);
    
    % Close the file.
    close(vidObj);
    
function [ED,BW] = getEdges(I)
    %I=imgaussfilt(I,1);
    level = graythresh(I);  % Calculates the gray threshold --KLB
    BW = double(im2bw(I,level^0.1)); %BW = BW version of image
    BW = medfilt2(BW,[5 5]); % Remove background noise only the outline of the mice remains
    ED = edge(BW,'canny'); %Edge detection (works a little better than Sobel)
return;

