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
function RatHeadTrack_nose
    close all;
    %Genaralized photo read from directory (create a folder 'video' that stores all
    %videos)
    path = [pwd '/video/*.tif'];
    files = dir(path);


        % Prepare the new file.
        vidObj = VideoWriter('test.avi');
        open(vidObj);

        % Create an animation.
        set(gca,'nextplot','replacechildren');

        nose_coord = [];
        head_coord = [];

    for i = 1:length(files)
    % Runs through all of the pictures for a specific view, converts them to
    % black and white, edge detection, then converts them all into a video.     

        %if the string name contains the second video set break from for
        %loop -- UM
        if ~isempty(strfind(files(i).name,'_c002')) 
            count=count+1;
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
        ED = getEdges(I);

        %identify nose and head points --UM
         if i==1
             fig = figure(1); h = imshow(ED,[0,1]); hold on 
             [x_start,y_start] = getpts(fig); % hand pick nose and head coordinates
             %(double click ends selection process)
             plot(x_start,y_start,'g'), hold off;
             [nrows,ncols]=size(get(h,'cdata'));
             xdata=get(h,'xdata'); ydata=get(h,'ydata');
             %Convert axes to pixel coordinates for accurate mapping. Was
             %not essential to the given data set
             x_nose = axes2pix(ncols,xdata,x_start(1)); y_nose = axes2pix(nrows,ydata,y_start(1));
             x_head = axes2pix(ncols,xdata,x_start(2)); y_head=axes2pix(nrows,ydata,y_start(2));
             nose=[x_nose y_nose]; head=[x_head,y_head];
             nose_coord = [nose_coord;nose]; head_coord=[head_coord;head];
         end

        %find the edges in the next image
        D_ED=getEdges(I_next);
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
        
        [~,ind_nose]=min(d_nose); [~,ind_head]=min(d_head);
        x_nose = E(ind_nose,1); y_nose = E(ind_nose,2);
        nose=[x_nose y_nose]; nose_coord =[nose_coord;nose];
        x_head = E(ind_head,1); y_head=E(ind_head,2);
        head=[x_head y_head]; head_coord=[head_coord;head];
        
        %Display nose and head points on image
        figure(3); imshow(ED,[0,1]), hold on;
        plot(nose_coord(i,1), nose_coord(i,2),'r*',...
            head_coord(i,1),head_coord(i,2),'go');
        hold off

        frame = getframe;
        writeVideo(vidObj,frame);
    end

    frame = getframe;
    writeVideo(vidObj,frame);
    
    % Close the file.
    close(vidObj);
    
function ED = getEdges(I)
    level = graythresh(I);  % Calculates the gray threshold --KLB
    BW = double(im2bw(I,level^0.5)); %BW = BW version of image
    BW = medfilt2(BW,[5 5]); % Remove background noise only the outline of the mice remains
    ED = edge(BW,'canny'); %Edge detection (works a little better than Sobel)
return;

