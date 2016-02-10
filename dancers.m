% Initialises general variables
Frames = 210;
Ignore_pos = [121,222; 150,223];
Threshold = 12;
Bins = 10;
Sep = [0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8];
Bin_edges = {Sep Sep };
Colour_thresh = 0.14;

Dancers = 4;
Colours = cellstr(['r.'; 'b.'; 'g.'; 'c.']);

all_observations = zeros(Dancers, Frames, 2);

% Initialises arrays containing the properties of each tracker (dancer)
avg_colour = zeros(Dancers, Bins, Bins);
avg_size = zeros(Dancers, 1);
matched = zeros(Dancers);


% First frame server as initialisation
Imstart = imread(strcat('DATA1/frame110.jpg'),'jpg');
% Constructs normalised red and green colour channels
redChannel = double(Imstart(:, :, 1));
greenChannel = double(Imstart(:, :, 2));
blueChannel = double(Imstart(:, :, 3));
rgb_sum = redChannel + greenChannel + blueChannel;
rgb_sum(rgb_sum == 0) = 1;
red = redChannel ./ rgb_sum;
green = greenChannel ./ rgb_sum;

[mask, candidate_regions] = getRegion( double(Imstart), Threshold, 0);
next_dancer = 1;

figure(1)
imshow(Imstart)
hold on
% Considers regions ordered by their size (starting with the largest) until
% Sufficiently many regions have been encountered
for region = 1 : size(candidate_regions)
    % Only considers regions that are not too close to excluded points
    if ~any(pdist2(Ignore_pos, [candidate_regions(region).Centroid(1), candidate_regions(region).Centroid(2)]) < 50)
        
        % Initialises size of the dancer
        avg_size(next_dancer) = candidate_regions(region).Area;
        
        % Initialises colour histogram of the dancer, by applying a mask to
        % the normalised red and green colourchannels
        col = (mask == region);        
        red_local = red(col);
        green_local = green(col);
        avg_colour(next_dancer,:,:) = hist3([red_local, green_local], 'Edges', Bin_edges) / avg_size(next_dancer);
        
        % Initialises next tracker here
        all_observations(next_dancer, 1, :) = [floor(candidate_regions(region).Centroid(1)), floor(candidate_regions(region).Centroid(2))];
        
        % Plots position of the tracker
        radius = sqrt(candidate_regions(region).Area / pi);
        for c = -0.99*radius: radius/10 : 0.99*radius
            r = sqrt(radius^2-c^2);
            plot(all_observations(next_dancer, 1,1) + c, all_observations(next_dancer, 1,2) + r, char(Colours(next_dancer)))
            plot(all_observations(next_dancer, 1,1) + c, all_observations(next_dancer, 1,2) - r, char(Colours(next_dancer)))
        end
        next_dancer = next_dancer + 1;
    end
    % Once all trackers have been assigned stop
    if next_dancer > Dancers
        break
    end
end
pause(1)


% Main program iterating over all frames and computing the trackers
% positions
for t = 2 : Frames
    disp(t)
    % Preparation - compute normalised red and green channels and extract region centres
    Imwork = imread(strcat('DATA1/frame', int2str(t + 109), '.jpg'),'jpg');
    redChannel = double(Imwork(:, :, 1));
    greenChannel = double(Imwork(:, :, 2));
    blueChannel = double(Imwork(:, :, 3));

    rgb_sum = redChannel + greenChannel + blueChannel;
    rgb_sum(rgb_sum == 0) = 1;

    red = redChannel ./ rgb_sum;
    green = greenChannel ./ rgb_sum;
    
    [~, properties] = getRegion( double(Imwork), Threshold, 0);
    candidate_regions = struct2cell(properties);
    centres = cell2mat(candidate_regions(2,:)');


    % for each tracker maintain a list of regiosn it can be matched with
    possible_regions = cell(Dancers, 1);
    % start by searching the immediate surroundings of the old position
    Search_Radius = 5;
    % as long as we haven't been able to find a suitable match keep
    % increasing the search readius until time-out
    attempts = 0;
    matched = 0;
    assignment = zeros(Dancers,1);
    while (~matched) && (attempts < 15)
        Search_Radius = Search_Radius + 10;
        % For every tracker (dancer) find regions within search radius of
        % old position
        for tracker = 1:Dancers
            distances = pdist2(centres, [all_observations(tracker, t-1,1), all_observations(tracker, t-1,2)]);
            % add newly encountered regions
            possible_regions(tracker,:) = {union(find(distances < Search_Radius)', possible_regions{tracker,:})};
        end

        % Compute all possible matches between regions and trackers as the
        % cartesian product of possible assignments for the individual
        % trackers
        [a,b,c,d] = ndgrid(possible_regions{:});
        cartProd = [a(:) b(:) c(:) d(:)];
        for m = 1 : size(cartProd, 1)
            % Reset the mask for the currently considered assignment
            [mask, props] = getRegion( double(Imwork), Threshold, 0);
            
            % Checks if regions need to be subdivided and performs a 
            % colour test on the resulting regions


            % Check if a single region is chosen multiple times and alter mask
            regions = unique(cartProd(m,:));
            % Variables for testing the final matchups against averag size
            % and colour of the respective trackers
            color_pass = 2;
            size_pass = 2;
            
            % If the trackers weren't matched to 4 different regions we
            % need to split at least one of them
            if ~(size([regions],2) == Dancers)
                % for every tracker, check if they have a region in common
                % with another tracker
                for track = 1 : Dancers
                    % If the tracker does share its assigned region
                    if sum(cartProd(m,:) == cartProd(m,track)) > 1
                        to_split = cartProd(m,track); % label of region to be split
                        [tmp_mask, new_min, new_max] = split_region(mask, to_split, sum(cartProd(m,:) == cartProd(m,track))); % resulting temporary mask and new region labels of even split
                        
                        split_participants = find(cartProd(m,:) == to_split); % tracker/dancer numbers assigned to the splitted region
                        allocations = perms(split_participants); % possible allocations to the resulting subregions
                        best_perm = zeros(numel(split_participants),1); % the best assignment to sub-regions
                        best_value = inf; % Best overall colour difference
                        
                        % for every permutation
                        for perm = 1 : size(allocations, 1)
                            val = 0;
                            
                            % Go through the assignments between sub-regions and their assigned tracker
                            % add the colour differences to 'val' 
                            for reg = 1 : size(allocations, 2)
                                col = (tmp_mask == reg + new_min - 1);
                                local_red = red(col);
                                local_green = green(col);
                                val = val + compareHists(hist3([local_red, local_green], 'Edges', Bin_edges) / sum(sum(col)), squeeze(avg_colour(allocations(perm, reg),:,:)));
                            end
                            % If it is better than the previous optimum
                            % remember it
                            if val < best_value
                                best_perm = allocations(perm, :);
                                best_value = val;
                            end
                        end
                        % Reassign the dancers according to the best
                        % permutation found and weight the final split by
                        % their average size
                        weights = zeros(size(split_participants,1),1);
                        for d = 1 : Dancers
                            if cartProd(m,d) == to_split
                                cartProd(m,d) = new_min - 1 + find(best_perm == d);
                                weights(find(best_perm == d)) = avg_size(d);
                            end
                        end
                        weights = weights / sum(weights);
                        [mask, ~, ~] = split_region_weighted(mask, to_split, weights); % final mask and new region labels
                    end
                    % Compute the colour difference between the dancer and
                    % its assigned region and test that it is within the
                    % allowed tolerance
                    col = (mask == cartProd(m,track));
                    local_red = red(col);
                    local_green = green(col);
                    col_shift = compareHists(hist3([local_red, local_green],'Edges', Bin_edges) / sum(sum(col)), squeeze(avg_colour(track,:,:)));
                    color_pass = color_pass - (col_shift > Colour_thresh);
                    % Do the same for the regions size
                    region_size = sum(sum(col));
                    size_pass = size_pass - (region_size > avg_size(track) * 1.3 || region_size < avg_size(track) * 0.7);                    

                end
            else
                % If all trackers were assigned to a different region
                % simply test the assignment with regards to region colours
                % and sizes
                for r = 1 : Dancers
                    % Colour comparison
                    col = (mask == cartProd(m,r));
                    local_red = red(col);
                    local_green = green(col);
                    col_shift = compareHists(hist3([local_red, local_green],'Edges', Bin_edges) / sum(sum(col)), squeeze(avg_colour(r,:,:)));
                    color_pass = color_pass - (col_shift > Colour_thresh);
                    % Size comparison
                    region_size = sum(sum(col));
                    size_pass = size_pass - (region_size > avg_size(r) * 1.3 || region_size < avg_size(r) * 0.7);      
                end
                
                % If all tests were successful update the tracker colours
                % and sizes
                if (color_pass == 2) && (size_pass == 2)
                    for r = 1 : Dancers
                        col = (mask == cartProd(m,r));
                        local_red = red(col);
                        local_green = green(col);
                        avg_colour(r,:,:) = hist3([local_red, local_green],'Edges', Bin_edges) / sum(sum(col)) / 5 + squeeze(avg_colour(r,:,:)) * 4 / 5 ;
                        avg_size(r) = sum(sum(col)) / 5 + avg_size(r) * 4 / 5;
                    end
                end
            end
            % If the colour test failed
            if color_pass < 1
                continue
            end
            % If the colour test was succesful
            if size_pass < 1
                continue
            end
            % update the assignment and leave
            assignment(:) = cartProd(m,:);
            matched = 1;
            break
        end
        attempts = attempts + 1;
    end
    
    
    if ~matched
        [mask, props] = getRegion( double(Imwork), Threshold, 0);

        for r = 1 : Dancers
            distances = pdist2(centres, [all_observations(r, t-1,1), all_observations(r, t-1,2)]);
            [min_Dist, minInd] = min(distances);
            assignment(r) = mask(floor(centres(minInd,2)), floor(centres(minInd,1)));
        end
        regions = unique(assignment);

        % If the trackers weren't matched to 4 different regions we
        % need to split at least one of them
        if ~(size([regions],2) == Dancers)
            % for every tracker, check if they have a region in common
            % with another tracker
            for track = 1 : Dancers
                % If the tracker does share its assigned region
                if sum(assignment == assignment(track)) > 1
                    to_split = assignment(track); % label of region to be split
                    [mask, new_min, new_max] = split_region(mask, to_split, sum(assignment == assignment(track))); % resulting mask and new region labels
                        
                    split_participants = find(assignment == to_split); % tracker/dancer numbers assigned to the splitted region
                    
                    allocations = perms(split_participants); % possible allocations to the resulting subregions
                    best_perm = zeros(numel(split_participants),1); % the best assignment to sub-regions
                    best_value = inf; % Best overall colour difference
                        
                    % for evey permutation
                    for perm = 1 : size(allocations, 1)
                        val = 0;

                        % Go through the assignments between sub-regions and their assigned tracker
                        % add the colour differences to 'val' 
                        for reg = 1 : size(allocations, 2)
                            col = (mask == reg + new_min - 1);
                            local_red = red(col);
                            local_green = green(col);
                            val = val + compareHists(hist3([local_red, local_green], 'Edges', Bin_edges) / sum(sum(col)), squeeze(avg_colour(allocations(perm, reg),:,:)));
                        end
                            % If it is better than the previous optimum
                            % remember it
                        if val < best_value
                            best_perm = allocations(perm, :);
                            best_value = val;
                        end
                    end
                    % Reassign the dancers according to the best
                    % permutation found
                    for d = 1 : Dancers
                        if assignment(d) == to_split
                            assignment(d) = new_min - 1 + find(best_perm == d);
                        end
                    end
                end
            end
        end
    end
    % Finalise the assignment
    hold off
    imshow(Imwork)
    hold on
    for tracker = 1 : Dancers
        [row,col] = find(mask == assignment(tracker));
        centre = [mean(col), mean(row)];
        all_observations(tracker, t, :) = centre;
        radius = sqrt(size(row,1) / pi);
        for c = -0.99*radius: radius/10 : 0.99*radius
            r = sqrt(radius^2-c^2);
            plot(all_observations(tracker, t,1) + c, all_observations(tracker, t,2) + r, char(Colours(tracker)))
            plot(all_observations(tracker, t,1) + c, all_observations(tracker, t,2) - r, char(Colours(tracker)))
        end
    end 
    pause(1)
end


