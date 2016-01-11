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
    BW1 = double(im2bw(I,level)); % BW1 = BW1 version of image (with lower grey level)
    BW2 = double(im2bw(I,level^0.5)); % BW2 = BW2 version of imag (with higher gery level)
    ED1 = edge(BW1,'Sobel'); % Edge detection 1
    ED2 = edge(BW2,'Sobel'); % Edge detection 2
    % imshow(ED);
    imshowpair(ED1,ED2,'montage') % Personally I think the BW2 is more suitable for detection -- CG

    frame = getframe;
    writeVideo(vidObj,frame);
end
    
    % Close the file.
    close(vidObj);
