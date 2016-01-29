% set(gca,'nextplot','replacechildren');
% Rikki Irwin

for i = 1:3
    if i>=1 && i <=9,
        str = '_c002s000100000';
    elseif i>=10 && i<= 99,
        str = '_c002s00010000';
    else
        str = '_c002s0001000';
    end

    I=imread(['video\',str,num2str(i),'.tif']);
    % Canny edge detection
    bw = edge(I,'canny',.8);
    
    % Hough transform, find peaks, and turn into lines. All the paramters
    % can be adjusted for how long a line needs to be, fillGap for merging
    % adjacent lines, how many peaks to look for, etc.
    [H,theta,rho] = hough(bw);
    P = houghpeaks(H,20,'threshold',ceil(.1*max(H(:))));
    lines = houghlines(bw,theta,rho,P,'FillGap',40,'MinLength',85);
    
    hold off
    imshow(I)
    hold on
    max_len = 0;
    
% Plot all the Hough lines    
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');


   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end

% Plot the mean of the two hough lines
xCen = 0;
yCen = 0;
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    cen = ceil([mean(xy(:,1)), mean(xy(:,2))]);
    xCen = xCen + cen(1);
    yCen = yCen + cen(2);
end
plot(xCen/k,yCen/k,'x','LineWidth',2,'Color','blue')
headCoords = [xCen/k, yCen/k];

% Find and plot the nose coordinates by taking the average of the points on
% the lines that are closest to each other
minDist = 1000000000;
p1 = [lines(1).point1; lines(1).point2];
p2 = [lines(2).point1; lines(2).point2];
for k = 1:2
    for j = 1:2
    dist = (p1(k,1) - p2(j,1))^2 + (p1(k,2)-p2(j,2))^2;
    if (dist < minDist)
        minDist = dist;
        ind1 = k;
        ind2 = j;
    end
    end
end

noseCoords = mean([p1(ind1,:); p2(ind2,:)]);
plot(noseCoords(1),noseCoords(2),'x','LineWidth',2,'Color','blue')
plot([noseCoords(1) headCoords(1)],[noseCoords(2) headCoords(2)],'LineWidth',2,'Color','blue')

% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');

pause(.01)

end
