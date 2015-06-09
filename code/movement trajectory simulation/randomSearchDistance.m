function totalDistance = randomSearchDistance(start,target,detectionRadius,movementType,mu,l_min,showPlot)
%   randomSearchDistance(start,target,detectionRadius,movementType,mu)
%   simulations and returns the distance traversed in the search for a 
%   target in a square with side length 1 (left bottom corner located [0 0]
%   The inputs are:
%        start               Start location, two value matrix: [x y]
%        target              Target location, two value matrix: [x y]
%        detectionRadius     Value > 0 indicating detection radius
%        movementType        Possible values:
%                               ballistic, levy, spiral, square, zigzag
%        mu                  Used for Levy walks. mu >= 3 is Brownian
%        l_min               Indicates the minimum step length (e.g., 0.05)
%        showPlot            Possible values: 0/false or 1/true
%                            0/false: no plot is generated
%                            1/true: a plot of the trajectory is generated


switch movementType
    case 'ballistic'
        currentLocation = start;
    case 'levy'
        currentLocation = start;
    case 'spiral'
        currentLocation = [.5 .5];
    case 'square'
        currentLocation = [.5 .5];
    case 'parallel'
        a = sqrt(detectionRadius ^2 / 2);
        currentLocation = [a a];
    case 'zigzag'
        currentLocation = [0 0];
    otherwise
        error('Unknown movement type.');
end


trajectory = currentLocation;

found = 0;
while(~found)
    
    if strcmp(movementType,'square') && size(trajectory,1) == 1
        currentLocation = [.5 .5];
        trajectory = currentLocation;
        step_lengths = [];
        directions = [];
    end
                
    if (sqrt(sum((currentLocation - target).^2)) - detectionRadius) <= eps
        found = 1;
    else
        
        switch movementType
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'ballistic' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                x = currentLocation(1);
                y = currentLocation(2);
                %If on a border, there are some angles that cannot be
                %chosen
                if x ~= 0 && x~= 1 && y ~= 0 && y~= 1
                    theta = rand*360;
                else
                    angleFound = 0;
                    while ~angleFound
                        theta = rand*360;
                        if ~((x <= 0 && theta >= 90 && theta <= 270) || ... %x = 0, no 90 - 270
                                (x >= 1 && ((theta >= 270 && theta <= 360) || (theta >= 0 && theta <= 90))) || ... %x = 1, no 270 - 90
                                (y <= 0 && ((theta >= 180 && theta <= 360) || (theta == 0))) || ... %y = 0, no 180 - 360
                                (y >= 1 && ((theta >= 0 && theta <= 180) || (theta == 360)))) %y = 1, no 0 - 180
                            angleFound = 1;
                        end
                    end
                end
                    
                %y = mx+b
                m = tand(theta);
                b = y - x*m;
                
                possible_next_coordinates = [Inf Inf; Inf Inf; Inf Inf; Inf Inf];
                %%Collision with bottom
                if y == 0
                    distance_bottom = Inf;
                else
                    y_bottom = 0;
                    x_bottom = -b/m;
                    possible_next_coordinates(1,1:2) = [x_bottom, y_bottom];
                    if x_bottom > 1 || x_bottom < 0
                        distance_bottom = Inf;
                    else
                        distance_bottom = sqrt((x-x_bottom)^2 + (y - y_bottom)^2);
                    end
                end
                    
                %%Collision with top
                if y == 1
                    distance_top = Inf;
                else
                    y_top = 1;
                    x_top = (1-b)/m;
                    possible_next_coordinates(2,1:2) = [x_top, y_top];
                    if x_top > 1 || x_top < 0
                        distance_top = Inf;
                    else
                        distance_top = sqrt((x-x_top)^2 + (y-y_top)^2);
                    end
                end
                    
                    
                %%Collision with left
                if x == 0
                    distance_left = Inf;
                else
                    y_left = b;
                    x_left = 0;
                    possible_next_coordinates(3,1:2) = [x_left, y_left];
                    if y_left > 1 || y_left < 0
                        distance_left = Inf;
                    else
                        distance_left = sqrt((x-x_left)^2 + (y-y_left)^2);
                    end
                end
                
                %Collision with right
                if x == 1
                    distance_right = Inf;
                else
                    y_right = m + b;
                    x_right = 1;
                    possible_next_coordinates(4,1:2) = [x_right, y_right];
                    if y_right > 1 || y_right < 0
                        distance_right = Inf;
                    else
                        distance_right = sqrt((x-x_right)^2 + (y-y_right)^2);
                    end
                end
                
                
                distances = [distance_bottom distance_top distance_left distance_right];
                index = find(distances==min(distances));
                nextLocation = possible_next_coordinates(index(1),:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'levy' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                theta = NaN;
                locationFound = 0;
                while ~locationFound
                    theta = rand*360;
%                     stepLength = l_min + (C*rand).^(1/(1-mu));
                    stepLength = l_min * rand.^(1/(1-mu));
                    nextLocation = currentLocation + [sind(theta)*stepLength cosd(theta)*stepLength];
                    locationFound = inpolygon(nextLocation(1),nextLocation(2), [0 1 1 0 0], [0 0 1 1 0]);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'spiral' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                a = 0;
                b = detectionRadius/(pi);

                theta = linspace(0,(1/detectionRadius)*2*pi,(1/detectionRadius)*10000);
                r = a + b*theta;
                xs = .5 + (r.*cos(theta))';
                ys = .5 + (r.*sin(theta))';
                locations = [xs ys];
                distances = sqrt(sum((locations - repmat(target,size(xs,1),1)).^2,2));
                indicesInRange = find(distances <= detectionRadius);
                trajectory = locations(1:indicesInRange(1),:);
                found = 1;
            case 'zigzag' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if currentLocation(1) == 0
                    nextLocation = [1 currentLocation(2) + detectionRadius];
                else
                    nextLocation = [0 currentLocation(2) + detectionRadius];
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'parallel'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                a = sqrt(detectionRadius ^2 / 2);
                if size(trajectory,1)==1 || (trajectory(end - 1,1) == a && currentLocation(1) == a)
                    nextLocation = [1-a currentLocation(2)];
                else if size(trajectory,1)==2 || (currentLocation(1) == 1-a && trajectory(end-1,1) == a)
                        nextLocation = [1-a currentLocation(2) + 2*a];
                    else if currentLocation(1) == 1-a && trajectory(end-1,1) == 1-a
                            nextLocation = [a currentLocation(2)];
                        else
                            nextLocation = [a currentLocation(2) + 2*a];
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'square'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                s = sqrt(2)*detectionRadius;
                stepLength = s*ceil(size(trajectory,1)/2);
                switch mod(size(trajectory,1),4)
                    case 0
                        nextLocation = [currentLocation(1) - stepLength currentLocation(2)];
                    case 1
                        nextLocation = [currentLocation(1) currentLocation(2)+stepLength];
                    case 2
                        nextLocation = [currentLocation(1) + stepLength currentLocation(2)];
                    case 3
                        nextLocation = [currentLocation(1) currentLocation(2) - stepLength];
                    otherwise
                        error('Something wrong with square.');
                end
                
                
            otherwise
                error('Unknown movement type.');
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%Test for encountering along path from currentLocation to
        %%nextLocation
        if ~strcmp(movementType,'spiral')
            x1 = currentLocation(1);
            x2 = nextLocation(1);
            y1 = currentLocation(2);
            y2 = nextLocation(2);
            x0 = target(1);
            y0 = target(2);
            distance_start_to_end = sqrt((x2-x1)^2 + (y2-y1)^2);
            distance_target_to_line = abs( (y2-y1)*x0 - (x2-x1)*y0 -x1*y2 + x2*y1 ) / distance_start_to_end;
            if distance_target_to_line <= detectionRadius
                %calculate if other distance is between 0 and positive the
                %distance between E and S
                distance_start_to_target = sqrt((x0-x1)^2 + (y0-y1)^2);
                
                theta_temp = asind(distance_target_to_line / distance_start_to_target);
                check = sum((target - currentLocation) .* (nextLocation - currentLocation));
                
                distance_pf = sqrt(detectionRadius^2 - distance_target_to_line^2);
                distance_ps = sqrt(distance_start_to_target^2 - distance_target_to_line^2);
                distance_sf = distance_ps - distance_pf;
                if check >= 0 && ...
                        distance_sf <= distance_start_to_end
                    nextLocation = currentLocation + distance_sf * (nextLocation - currentLocation) / distance_start_to_end;
                    found = 1;
                end
            end
            
            currentLocation = nextLocation;
            trajectory = [trajectory ; nextLocation];
        end
    end
    
    if found
        currentLocation = target;
        trajectory = [trajectory ; target];
    end
end

legDistances = sqrt(sum((trajectory(2:end,:) - trajectory(1:(end-1),:)).^2,2));
totalDistance = sum(legDistances);

if showPlot
    close all
    figure
    plot([0 1 1 0 0],[0 0 1 1 0],'r-');
    hold on
    plot(trajectory(:,1),trajectory(:,2),'k-');
    xlim([-.1,1.1]);
    ylim([-.1,1.1]);
    plot(target(1),target(2),'r*');
    plot(currentLocation(1,1),currentLocation(1,2),'go');
    title(['Number of steps: ' num2str(size(trajectory,1)-1) ', dection radius = ' num2str(detectionRadius)]);
end
