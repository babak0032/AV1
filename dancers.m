Dancers = 4;
Colours = cellstr(['r.'; 'b.'; 'g.'; 'c.']);
Frames = 210;
Particles = 2;
Ignore_pos = [121,222; 150,223];
Threshold = 12;
Bins = 10;
Sep = [0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8];
Bin_edges = {Sep Sep };%{[0, 5, 15, 30, 50, 75] [0, 5, 15, 30, 50, 75]};
Colour_thresh = 0.14; % 0.1554;


Imback = double(imread('DATA1/bgframe.jpg','jpg'));
[MR,MC,Dim] = size(Imback);
% Tracker for each Dancer
avg_colour = zeros(Dancers, Bins, Bins);
avg_size = zeros(Dancers, 1);
matched = zeros(Dancers);

all_observations = zeros(Dancers, Frames, 2);
path = zeros(Dancers, Frames, 2);
track_pos = zeros(Dancers, Particles, Frames, 4);
track_weights = zeros(Dancers, Particles);
track_P = zeros(Dancers, Particles, 4, 4);

track_P(:,:,1,1) = 5;
track_P(:,:,2,2) = 5;
track_P(:,:,3,3) = 0;
track_P(:,:,4,4) = 0;

Imstart = imread(strcat('DATA1/frame110.jpg'),'jpg');
next_dancer = 1;
[mask, candidate_regions] = getRegion( double(Imstart), Threshold, 0);
figure(1)
imshow(Imstart)

for region = 1 : size(candidate_regions)
    if ~any(pdist2(Ignore_pos, [candidate_regions(region).Centroid(1), candidate_regions(region).Centroid(2)]) < 50)

        avg_size(next_dancer) = candidate_regions(region).Area;

        redChannel = double(Imstart(:, :, 1));
        greenChannel = double(Imstart(:, :, 2));
        blueChannel = double(Imstart(:, :, 3));

        rgb_sum = redChannel + greenChannel + blueChannel;
        rgb_sum(rgb_sum == 0) = 1;

        red = redChannel ./ rgb_sum;
        green = greenChannel ./ rgb_sum;
        
        col = (mask == region);        
        red = red(col);
        green = green(col);
        avg_colour(next_dancer,:,:) = hist3([red, green], 'Edges', Bin_edges) / avg_size(next_dancer);
        
        for k = 1 : Particles
              track_pos(next_dancer, k,1,:) = [floor(candidate_regions(region).Centroid(1) + 15*rand(1)-7), floor( candidate_regions(region).Centroid(2) + 15*rand(1)-7),0,0];
              track_weights(next_dancer, k) = 1 / Particles;
        end
        path(next_dancer, 1, :) = [floor(candidate_regions(region).Centroid(1)), floor(candidate_regions(region).Centroid(2))];
        all_observations(next_dancer, 1, :) = [floor(candidate_regions(region).Centroid(1)), floor(candidate_regions(region).Centroid(2))];
        radius = sqrt(candidate_regions(region).Area / pi);
        next_dancer = next_dancer + 1;

    end
    if next_dancer > Dancers
        break
    end
end
hold off
figure(1)



for t = 2 : Frames
    disp(t)
    % Preparation - extract region centres
    time_step = 0;
    Search_Radius = 5;
    Imwork = imread(strcat('DATA1/frame', int2str(t + 109), '.jpg'),'jpg');
    
    [~, properties] = getRegion( double(Imwork), Threshold, 0);
    candidate_regions = struct2cell(properties);
    centres = cell2mat(candidate_regions(2,:)');
    matched = 0;

    % for each tracker find regions close to estimated position
    possible_regions = cell(Dancers, 1);
    
    % as long as we haven't been able to find a suitable match keep
    % increasing the time step
    while ~all(matched)
        time_step = time_step + 1;
        Search_Radius = Search_Radius + 10;
        hold off
        hold on
        imshow(Imwork)       

        
        estimates = zeros(Dancers, 2);
        
        % For every tracker (dancer) predict new positions and find the
        % regions close to it
        for tracker = 1:Dancers
            [new_state, temp_state] = predict(squeeze(track_pos(tracker, :, t - 1, : )), squeeze(track_P(tracker, :, :, :)), track_weights(tracker, : ), time_step);
            estimates(tracker, :) = [new_state(1), new_state(2)];
            distances = pdist2(centres, estimates(tracker, :));
            possible_regions(tracker,:) = {union(find(distances < Search_Radius)', possible_regions{tracker,:})};
            
%             for k = 1 : Particles
%               plot(temp_state(k,1), temp_state(k,2), char(Colours(tracker)))
%             end
            
        end

        for tracker = 1:Dancers
            radius = Search_Radius;
            for c = -0.99*radius: radius/10 : 0.99*radius
              r = sqrt(radius^2-c^2);
              plot(estimates(tracker,1) + c, estimates(tracker,2) + r, char(Colours(tracker)))
              plot(estimates(tracker,1) + c, estimates(tracker,2) - r, char(Colours(tracker)))
            end
        end
%         for tracker = 1:Dancers
%             radius = Search_Radius;
%             for c = -0.99*radius: radius/10 : 0.99*radius
%               r = sqrt(radius^2-c^2);
%               plot(path(tracker,t-1, 1) + c, path(tracker, t-1, 2) + r, char(Colours(tracker)))
%               plot(path(tracker,t-1, 1) + c, path(tracker, t-1, 2) - r, char(Colours(tracker)))
%             end
%         end
        pause(0.05)
        % Compute all possible matches between regions and trackers
        [a,b,c,d] = ndgrid(possible_regions{:});
        cartProd = [a(:) b(:) c(:) d(:)];
        for m = 1 : size(cartProd, 1)

            %alsdfgh = cartProd(m,:)
            [mask, props] = getRegion( double(Imwork), Threshold, 0);
            
            % Checks if regions need to be subdivided and performs a simple
            % colour test on the regions
            redChannel = double(Imwork(:, :, 1));
            greenChannel = double(Imwork(:, :, 2));
            blueChannel = double(Imwork(:, :, 3));
            
            rgb_sum = redChannel + greenChannel + blueChannel;
            rgb_sum(rgb_sum == 0) = 1;

            red = redChannel ./ rgb_sum;
            green = greenChannel ./ rgb_sum;

            % Check if region is occupied multiple times and alter mask
            regions = unique(cartProd(m,:));
            color_pass = 2;
            size_pass = 1;
            
            % If we don't have 4 different regions
            if ~(size([regions],2) == Dancers)
                disp('overcrowded')
                
                % for every tracker, check if they choose region in common
                % with another tracker
                for track = 1 : Dancers
                    
                    % If the tracker does choose a region in common
                    if sum(cartProd(m,:) == cartProd(m,track)) > 1
                        disp('Splitting')
                        % disp( cartProd(m, :) )
                        to_split = cartProd(m,track); % the region number to be split
                                           
                        [mask, new_min, new_max] = split_region(mask, to_split, sum(cartProd(m,:) == cartProd(m,track)));
                        
                        split_participants = find(cartProd(m,:) == to_split); % tracker/dancer numbers of the split
                        allocations = perms(split_participants); % permutations of the trackers who share the region
                        best_perm = zeros(numel(split_participants),1); % the best permutation (which tracker assigns to whcih sub-rgion)
                        best_value = inf; % Best colour difference sum (The lowest)
                        
                        % for evey premutation
                        for perm = 1 : size(allocations, 1)
                            val = 0;
                            
                            % Go through every sub-region and add its colour difference to 'val' 
                            for reg = 1 : size(allocations, 2)
                                col = (mask == reg + new_min - 1);
                                local_red = red(col);
                                local_green = green(col);
                                val = val + compareHists(hist3([local_red, local_green], 'Edges', Bin_edges) / sum(sum(col)), squeeze(avg_colour(allocations(perm, reg),:,:)));
                            end
                            if val < best_value
                                best_perm = allocations(perm, :);
                                best_value = val;
                            end
                        end
                        for d = 1 : Dancers
                            if cartProd(m,d) == to_split
                                cartProd(m,d) = new_min - 1 + find(best_perm == d);
                            end
                        end
                   
                    end
                    % disp('singleton')
                    col = (mask == cartProd(m,track));
                    local_red = red(col);
                    local_green = green(col);
                    col_shift = compareHists(hist3([local_red, local_green],'Edges', Bin_edges) / sum(sum(col)), squeeze(avg_colour(track,:,:)));
                    color_pass = color_pass - (col_shift > Colour_thresh);
                    % size_pass = size_pass & sum(sum(mask == cartProd(m,d))) < avg_size(track) * 1.2 & sum(sum(mask == cartProd(m,d))) > avg_size(track) * 0.8;                    

                end
            else
                
                % if size(cartProd,1) > 1
                    % Do simple colour test
                for r = 1 : Dancers
                    col = (mask == cartProd(m,r));
                    local_red = red(col);
                    local_green = green(col);
                    col_shift = compareHists(hist3([local_red, local_green],'Edges', Bin_edges) / sum(sum(col)), squeeze(avg_colour(r,:,:)));
                    color_pass = color_pass - (col_shift > Colour_thresh);
                    % size_pass = size_pass & sum(sum(mask == cartProd(m,r))) < avg_size(r) * 1.2 & sum(sum(mask == cartProd(m,r))) > avg_size(r) * 0.8;
                end
                % end

                if color_pass == 2
                    for r = 1 : Dancers
                        % TODO update
                        col = (mask == cartProd(m,r));
                        local_red = red(col);
                        local_green = green(col);
                        avg_colour(r,:,:) = hist3([local_red, local_green],'Edges', Bin_edges) / sum(sum(col)) / 5 + squeeze(avg_colour(r,:,:)) * 4 / 5 ;
                        % avg_size(r) = sum(sum(mask == cartProd(m,r)));
                    end
                end
                
                
            end
            
            % this just calculates the colour difference (For debugging
            % purposes)
            values = unique(mask)';
            col_diffs = zeros(1,Dancers);% zeros(Dancers, size(values,2) -1 );
            for d = 1 : Dancers
                % for reg = 1 : size(values,2) -1
                col = (mask == cartProd(m,d)); % values(reg + 1));
                local_red = red(col);
                local_green = green(col);
                col_diffs(1,d) = compareHists(hist3([local_red, local_green],'Edges', Bin_edges) / sum(sum(col)), squeeze(avg_colour(d,:,:)));
                % end
            end
            
            col_diffs;
            
            col_diffs_2 = zeros(Dancers, size(values,2) -1 );
            for d = 1 : Dancers
                for reg = 1 : size(values,2) -1
                  col = (mask == values(reg + 1));
                  local_red = red(col);
                  local_green = green(col);
                  col_diffs_2(d,reg) = compareHists(hist3([local_red, local_green],'Edges', Bin_edges) / sum(sum(col)), squeeze(avg_colour(d,:,:)));
                end
            end           
            
            
            allocs = cartProd(m,:);
            col_diffs_2;
            pause(0.05)
            
            color_pass;
            if color_pass < 1
                imshow(label2rgb(mask))
                pause(1)
                disp('COLORS');
                continue
            end
%             if ~size_pass
%                 disp(avg_size)
%                 for d = 1 : Dancers
%                     sum(sum(mask == cartProd(m,d)))
%                 end
% 
%                 disp('SIZE');
%                 continue
%             end
            % Passed all test, update

            for r = 1 : Dancers
                [row,col] = find(mask == cartProd(m,r));
                centre = [mean(col), mean(row)];
                all_observations(r, t, :) = centre;
                [new_x, new_P, new_weights, pos] = condense_function(squeeze(track_pos(r,:,t-1,:)), squeeze(track_P(r, :, :, :)), squeeze(track_weights(r,:)), centre, time_step);
                track_pos(r, :, t, :) = new_x; 
                track_weights(r, :) = new_weights;
                track_P(r, :, :, :) = new_P;
                path(r, t, :) = [pos(1), pos(2)];
            end            
            matched = 1;
            break
        end
    end
end


