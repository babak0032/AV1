Dancers = 4;
Frames = 209;
Particles = 100;
Ignore_pos = [121,222; 150,223];

Imback = double(imread('DATA1/bgframe.jpg','jpg'));
[MR,MC,Dim] = size(Imback);
% Tracker for each Dancer
avg_colour = zeros(Dancers, 3, 25);
avg_size = zeros(Dancers,1);
matched = zeros(Dancers);

path = zeros(Dancers, Frames, 2);
track_pos = zeros(Dancers, Particles, Frames, 4); 
track_weights = zeros(Dancers, Particles);
track_P = zeros(Dancers, Particles, 4, 4);

track_P(:,:,1,1) = 5;
track_P(:,:,2,2) = 5;
track_P(:,:,3,3) = 5;
track_P(:,:,4,4) = 5;

% for i =1:Frames
%     [candidate_regions, mask] = getRegion( double(imread(strcat('DATA1/frame', int2str(i + 109), '.jpg'),'jpg')), 10, 0);
%     imshow(mask)
%     pause(0.05)
% end

Imstart = imread(strcat('DATA1/frame110.jpg'),'jpg');
next_dancer = 1;
mask = getRegion( double(Imstart), 10, 0);
candidate_regions = regions_sorted(mask);
figure(1)

imshow(Imstart)
hold on

disp(1)
for region = 1 : size(candidate_regions)
    if ~any(pdist2(Ignore_pos, [candidate_regions(region).Centroid(1), candidate_regions(region).Centroid(2)]) < 50)
        
        col = mask == region;
        
        redChannel = Imstart(:, :, 1);
        greenChannel = Imstart(:, :, 2);
        blueChannel = Imstart(:, :, 3);
        
        [redhist, ~] = imhist(redChannel(col), 25);
        [greenhist, ~] = imhist(greenChannel(col), 25);
        [bluehist, ~] = imhist(blueChannel(col), 25);
        
        avg_colour(next_dancer,:,:) = [redhist, greenhist, bluehist]';       
        avg_size(next_dancer) = candidate_regions(region).Area;
        
        for k = 1 : Particles
              track_pos(next_dancer, k,1,:) = [floor(candidate_regions(region).Centroid(1) + 15*rand(1)-7), floor( candidate_regions(region).Centroid(2) + 15*rand(1)-7),0,0];
              track_weights(next_dancer, k) = 1 / Particles;
        end
        path(next_dancer, 1, :) = [floor(candidate_regions(region).Centroid(1)), floor(candidate_regions(region).Centroid(2))];
        radius = sqrt(candidate_regions(region).Area / pi);
%         for k = 1 : Particles
%             plot(track_pos(next_dancer,k,1,1), track_pos(next_dancer, k,1,2),'r.')
%         end
        for c = -0.99*radius: radius/10 : 0.99*radius
          r = sqrt(radius^2-c^2);
          plot(candidate_regions(region).Centroid(1) + c, candidate_regions(region).Centroid(2) + r,'r.')
          plot(candidate_regions(region).Centroid(1) + c, candidate_regions(region).Centroid(2) - r,'r.')
        end
        next_dancer = next_dancer + 1;
        disp([candidate_regions(region).Centroid(1), candidate_regions(region).Centroid(2)])
    end
    if next_dancer > Dancers
        break
    end
end
hold off

pause(0.5)




for t = 2 : 32
    disp(t)
    % Preparation - extract region centres
    time_step = 0;
    Imwork = imread(strcat('DATA1/frame', int2str(t + 109), '.jpg'),'jpg');
    
    imshow(Imwork)
    hold on
    
    figure(2)
    mask = getRegion( double(Imwork), 10, 0);
    imshow(mask)
    figure(1)
    
    candidate_regions = struct2cell(regions_sorted(getRegion( double(Imwork), 10, 0)));
    centres = cell2mat(candidate_regions(2,:)');
    matched = 0;

    % as long as we haven't been able to find a suitable match keep
    % increasing the time step
    while ~all(matched)
        time_step = time_step + 1;
        % for each tracker find regions close to estimated position
        possible_regions = cell(Dancers, 1);
        
        estimates = zeros(Dancers, 2);
        
        for tracker = 1:Dancers
            [new_state, temp_state] = predict(squeeze(track_pos(tracker, :, t - 1, : )), track_weights(tracker,:), time_step);
            estimates(tracker, :) = [new_state(1), new_state(2)];
            distances = pdist2(centres, estimates(tracker, :));
            possible_regions(tracker,:) = {find(distances < 20)'};
            
            for k = 1 : Particles
              plot(temp_state(k,1), temp_state(k,2),'b.')
            end
            
            
        end

        for tracker = 1:Dancers
            radius = 30;
            for c = -0.99*radius: radius/10 : 0.99*radius
              r = sqrt(radius^2-c^2);
              plot(estimates(tracker,1) + c, estimates(tracker,2) + r,'g.')
              plot(estimates(tracker,1) + c, estimates(tracker,2) - r,'g.')
            end
        end
        for tracker = 1:Dancers
            radius = 30;
            for c = -0.99*radius: radius/10 : 0.99*radius
              r = sqrt(radius^2-c^2);
              plot(path(tracker,t-1, 1) + c, path(tracker, t-1, 2) + r,'r.')
              plot(path(tracker,t-1, 1) + c, path(tracker, t-1, 2) - r,'r.')
            end
        end
        pause(0.5)
        % Compute all possible matches between regions and trackers
        [a,b,c,d] = ndgrid(possible_regions{:});
        cartProd = [a(:) b(:) c(:) d(:)];
        for m = 1:size(cartProd, 1)
            
            mask = getRegion( double(Imwork), 10, 0);

%             % First test check if region is occupied twice and alter mask
%             regions = unique(cartProd(m,:));
%             if ~(size([regions],2) == Dancers)
%                 disp('Something')
%                 for track = 1 : size(cartProd(m,:))
%                     if sum(cartProd(m,:) == cartProd(m,track)) > 1
%                         disp('Splitting')
%                         points = estimates(cartProd(m,:) == cartProd(m,track), :);
%                         mask = split_region(mask, cartProd(m,track), points);
%                     end
%                 end
%             end
            
            
            % Second test - size
            %
            if ~all([candidate_regions{1, cartProd(m,:)}]' - 0.5 * avg_size(:) > 0 )
                continue
            end
            
            
            % Third test
            redChannel = Imwork(:, :, 1);
            greenChannel = Imwork(:, :, 2);
            blueChannel = Imwork(:, :, 3);
            color_pass = 1;
            
            for r = 1 : Dancers
                col = mask == cartProd(m,r);

                [redhist, ~] = imhist(redChannel(col), 25);
                [greenhist, ~] = imhist(greenChannel(col), 25);
                [bluehist, ~] = imhist(blueChannel(col), 25);

                col_shift = pdist2(redhist, squeeze(avg_colour(r,1,:))) + ...
                    pdist2(greenhist, squeeze(avg_colour(r,2,:))) + ...
                    pdist2(bluehist, squeeze(avg_colour(r,3,:)));
                color_pass = color_pass & (col_shift < 200);
                
            end
            if ~color_pass
                disp('COLORS')
                continue
            end
            for r = 1 : Dancers
                [new_x, new_P, new_weights, pos] = condense_function(squeeze(track_pos(r,:,t-1,:)), squeeze(track_P(r, :, :, :)), squeeze(track_weights(r,:)), centres(cartProd(m,r),:), time_step);
                track_pos(r, :, t, :) = new_x; 
                track_weights(r, :) = new_weights;
                track_P(r, :, :, :) = new_P;
                path(r, t, :) = [pos(1), pos(2)];
                disp([pos(1), pos(2)])
            end            
            matched = 1;
            break
        end
    end
end
        
        
        
        
        
%             if any()
%                 
%                 col_shift(distances > 100) = Inf;
%                 [~,choice] = min(col_shift);
%                 matched(tracker) = choice;
%                 %TODO update
%                 track_pos(tracker, 1,t,:) = [centres(choice,1), centres(choice,2),0,0,0,0];
%             else
%                 disp('FAIL')
%                 plot(new_state(1) + c, new_state(2) + r,'g.');
%             end
%         end
%     end
%     for tracker = 1:Dancers
%         radius = sqrt(candidate_regions{1, matched(tracker)} / pi);
% %         for k = 1 : Particles
% %             plot(track_pos(tracker,k,1,1), track_pos(tracker, k,1,2),'r.')
% %         end
%         for c = -0.99*radius: radius/10 : 0.99*radius
%           r = sqrt(radius^2-c^2);
%           plot(centres(matched(tracker),1) + c, centres(matched(tracker),2) + r,'r.')
%           plot(centres(matched(tracker),1) + c, centres(matched(tracker),2) - r,'r.')
%         end
%     end
%     hold off
%     pause(0.1)
% end
