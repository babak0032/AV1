Dancers = 4;
Colours = cellstr(['r.'; 'b.'; 'g.'; 'c.']);
Frames = 209;
Particles = 100;
Ignore_pos = [121,222; 150,223];

Imback = double(imread('DATA1/bgframe.jpg','jpg'));
[MR,MC,Dim] = size(Imback);
% Tracker for each Dancer
avg_colour = zeros(Dancers, 3, 256);

matched = zeros(Dancers);

path = zeros(Dancers, Frames, 2);
track_pos = zeros(Dancers, Particles, Frames, 4);
track_weights = zeros(Dancers, Particles);
track_P = zeros(Dancers, Particles, 4, 4);

track_P(:,:,1,1) = 5;
track_P(:,:,2,2) = 5;
track_P(:,:,3,3) = 0;
track_P(:,:,4,4) = 0;

% for i =1:Frames
%     [candidate_regions, mask] = getRegion( double(imread(strcat('DATA1/frame', int2str(i + 109), '.jpg'),'jpg')), 10, 0);
%     imshow(mask)
%     pause(0.05)
% end

Imstart = imread(strcat('DATA1/frame110.jpg'),'jpg');
next_dancer = 1;
[mask, candidate_regions] = getRegion( double(Imstart), 10, 0);
figure(2)
imshow(Imstart)

for region = 1 : size(candidate_regions)
    if ~any(pdist2(Ignore_pos, [candidate_regions(region).Centroid(1), candidate_regions(region).Centroid(2)]) < 50)
        
        col = mask == region;
        
        redChannel = Imstart(:, :, 1);
        greenChannel = Imstart(:, :, 2);
        blueChannel = Imstart(:, :, 3);
        
        [redhist, ~] = imhist(redChannel(col), 256);
        [greenhist, ~] = imhist(greenChannel(col), 256);
        [bluehist, ~] = imhist(blueChannel(col), 256);
        
        avg_colour(next_dancer,:,:) = [redhist, greenhist, bluehist]';
        
        for k = 1 : Particles
              track_pos(next_dancer, k,1,:) = [floor(candidate_regions(region).Centroid(1) + 15*rand(1)-7), floor( candidate_regions(region).Centroid(2) + 15*rand(1)-7),0,0];
              track_weights(next_dancer, k) = 1 / Particles;
        end
        path(next_dancer, 1, :) = [floor(candidate_regions(region).Centroid(1)), floor(candidate_regions(region).Centroid(2))];
        radius = sqrt(candidate_regions(region).Area / pi);
        
%         for c = -0.99*radius: radius/10 : 0.99*radius
%           r = sqrt(radius^2-c^2);
%           plot(candidate_regions(region).Centroid(1) + c, candidate_regions(region).Centroid(2) + r, 'r.')
%           plot(candidate_regions(region).Centroid(1) + c, candidate_regions(region).Centroid(2) - r, char(Colours(next_dancer)))
%         end
        next_dancer = next_dancer + 1;

    end
    if next_dancer > Dancers
        break
    end
end
pause(0.5)
hold off




for t = 2 : 32

    % Preparation - extract region centres
    time_step = 0;
    Imwork = imread(strcat('DATA1/frame', int2str(t + 109), '.jpg'),'jpg');
    
%     figure(1)
%     mask = getRegion( double(Imwork), 10, 1);
%     imshow(mask)
%     figure(2)
    [~, properties] = getRegion( double(Imwork), 10, 0);
    candidate_regions = struct2cell(properties);
    centres = cell2mat(candidate_regions(2,:)');
    matched = 0;

    % for each tracker find regions close to estimated position
    possible_regions = cell(Dancers, 1);
    
    % as long as we haven't been able to find a suitable match keep
    % increasing the time step
    while ~all(matched)
        time_step = time_step + 1;
        
        hold off
        imshow(Imwork)
        hold on

        
        estimates = zeros(Dancers, 2);
        
        for tracker = 1:Dancers
            [new_state, temp_state] = predict(squeeze(track_pos(tracker, :, t - 1, : )), squeeze(track_P(tracker, :, :, :)), track_weights(tracker, : ), time_step);
            estimates(tracker, :) = [new_state(1), new_state(2)];
            distances = pdist2(centres, estimates(tracker, :));
            possible_regions(tracker,:) = {union(find(distances < 60)', possible_regions{tracker,:})};
            
            for k = 1 : Particles
              plot(temp_state(k,1), temp_state(k,2), char(Colours(tracker)))
            end
            
        end

        for tracker = 1:Dancers
            radius = 60;
            for c = -0.99*radius: radius/10 : 0.99*radius
              r = sqrt(radius^2-c^2);
              plot(estimates(tracker,1) + c, estimates(tracker,2) + r, char(Colours(tracker)))
              plot(estimates(tracker,1) + c, estimates(tracker,2) - r, char(Colours(tracker)))
            end
        end
        for tracker = 1:Dancers
            radius = 60;
            for c = -0.99*radius: radius/10 : 0.99*radius
              r = sqrt(radius^2-c^2);
              plot(path(tracker,t-1, 1) + c, path(tracker, t-1, 2) + r, char(Colours(tracker)))
              plot(path(tracker,t-1, 1) + c, path(tracker, t-1, 2) - r, char(Colours(tracker)))
            end
        end
        % celldisp(possible_regions)
        if t > 18
            pause(0.1)
        else
            pause(0.1)
        end
        
        % Compute all possible matches between regions and trackers
        [a,b,c,d] = ndgrid(possible_regions{:});
        cartProd = [a(:) b(:) c(:) d(:)];
        for m = 1:size(cartProd, 1)
            %alsdfgh = cartProd(m,:)
            [mask, props] = getRegion( double(Imwork), 10, 0);
            

            % First test check if region is occupied twice and alter mask
            regions = unique(cartProd(m,:));
            if ~(size([regions],2) == Dancers)
                disp('Something')
                for track = 1 : size(cartProd(m,:))
                    if sum(cartProd(m,:) == cartProd(m,track)) > 1
                        disp('Splitting')
                        disp( cartProd(m,track) )
                        nr = size(cartProd(m,:) == cartProd(m,track));
                        [mask, index ] = split_region(mask, cartProd(m,track), props, 2);
                        for d = 1 : Dancers
                            if cartProd(m,d) == cartProd(m,track)
                                cartProd(m,d) = index;
                                index = index + 1;
                            end
                        end
                    end
                end
                imshow(label2rgb(mask))
                pause(2)
            end
            
            
            
            % Third test
            redChannel = Imwork(:, :, 1);
            greenChannel = Imwork(:, :, 2);
            blueChannel = Imwork(:, :, 3);
            color_pass = 1;
            
            col_diffs = zeros(Dancers, size(candidate_regions,2));
            for d = 1 : Dancers
                for reg = 1 : size(candidate_regions,2)
                    col = mask == reg;

                    [redhist, ~] = imhist(redChannel(col), 256);
                    [greenhist, ~] = imhist(greenChannel(col), 256);
                    [bluehist, ~] = imhist(blueChannel(col), 256);
                    
                    col_diffs(d,reg) = norm(redhist' - squeeze(avg_colour(d,1,:))') + ...
                        norm(greenhist' - squeeze(avg_colour(d,2,:))') + ...
                        norm(bluehist' - squeeze(avg_colour(d,3,:))');
                end
                [~,min_id] = min(col_diffs(d,:));
            end
            col_diffs;

            
            
            for r = 1 : Dancers
                col = mask == cartProd(m,r);

                [redhist, ~] = imhist(redChannel(col), 256);
                [greenhist, ~] = imhist(greenChannel(col), 256);
                [bluehist, ~] = imhist(blueChannel(col), 256);

                col_shift = norm(redhist' - squeeze(avg_colour(r,1,:))') + ...
                        norm(greenhist' - squeeze(avg_colour(r,2,:))') + ...
                        norm(bluehist' - squeeze(avg_colour(r,3,:))');
                col_shift;
                color_pass = color_pass & (col_shift < 300);
            end
            
            if ~color_pass
                disp('COLORS');
                continue
            end
            % Passed all test, update
            for r = 1 : Dancers
                [new_x, new_P, new_weights, pos] = condense_function(squeeze(track_pos(r,:,t-1,:)), squeeze(track_P(r, :, :, :)), squeeze(track_weights(r,:)), centres(cartProd(m,r),:), time_step);
                track_pos(r, :, t, :) = new_x; 
                track_weights(r, :) = new_weights;
                track_P(r, :, :, :) = new_P;
                path(r, t, :) = [pos(1), pos(2)];
                
                col = mask == cartProd(m,r);

                [redhist, ~] = imhist(redChannel(col), 256);
                [greenhist, ~] = imhist(greenChannel(col), 256);
                [bluehist, ~] = imhist(blueChannel(col), 256);
                
                %avg_colour(r,:,:) = [redhist, greenhist, bluehist]'/ t + squeeze(avg_colour(r,:,:)) * (t-1) / t ;
            end            
            matched = 1;
            break
        end
    end
end


