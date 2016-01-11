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


    % Prepare the new file.
    vidObj = VideoWriter('test.avi');
    open(vidObj);
 
    % Create an animation.
    set(gca,'nextplot','replacechildren');

for i = 1:140
% Runs through all of the pictures for a specific view, converts them to
% black and white, edge detection, then converts them all into a video.
    
    if i>=1 && i <=9,
        str = '_c001s000100000';
    elseif i>=10 && i<= 99,
        str = '_c001s00010000';
    else
        str = '_c001s0001000';
    end
    
    % Reads in the file   
    I=imread([str,num2str(i),'.tif']);
    
    
    level = graythresh(I);  % Calculates the gray threshold --KLB
    BW = double(im2bw(I,level)); %BW = BW version of image
    ED = edge(BW,'Sobel'); %Edge detection
    imshow(ED);
    frame = getframe;
    writeVideo(vidObj,frame);
end
    
    % Close the file.
    close(vidObj);