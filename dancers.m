Dancers = 4;
Frames = 209;
Particles = 100;
Ignore_pos = [121,222; 150,223];

Imback = double(imread('DATA1/bgframe.jpg','jpg'));
[MR,MC,Dim] = size(Imback);
% Tracker for each Dancer
avg_colour = zeros(Dancers, 3);
avg_size = zeros(Dancers);
matched = zeros(Dancers);

track_pos = zeros(Dancers, Particles, Frames, 6); 
track_weights = zeros(Dancers, Particles);
track_P = zeros(Dancers, Particles, 4, 4);

track_P(:,:,:,1,1) = 100;
track_P(:,:,:,2,2) = 100;
track_P(:,:,:,3,3) = 100;
track_P(:,:,:,4,4) = 100;

% for i =1:Frames
%     [candidate_regions, mask] = getRegion( double(imread(strcat('DATA1/frame', int2str(i + 109), '.jpg'),'jpg')), 10, 0);
%     imshow(mask)
%     pause(0.05)
% end
Imstart = imread(strcat('DATA1/frame110.jpg'),'jpg');
next_dancer = 1;
[candidate_regions, mask] = getRegion( double(Imstart), 10, 0);
figure(1)

imshow(Imstart)
hold on
for region = 1 : size(candidate_regions)
    if ~any(pdist2(Ignore_pos, [candidate_regions(region).Centroid(1), candidate_regions(region).Centroid(2)]) < 50)

        redChannel = int64(Imstart(:, :, 1));
        greenChannel = int64(Imstart(:, :, 2));
        blueChannel = int64(Imstart(:, :, 3));
        
        col = mask == region; % for example 2 or 3 or whatever region you want.
        
        avg_colour(next_dancer,:) = [mean(redChannel(col)), mean(greenChannel(col)), mean(blueChannel(col))];       
        avg_size(next_dancer) = candidate_regions(region).Area;
        
        for k = 1 : Particles
              track_pos(next_dancer, k,1,:) = [floor(candidate_regions(region).Centroid(1) + 15*rand(1)-7), floor( candidate_regions(region).Centroid(2) + 15*rand(1)-7),0,0,0,0];
              track_weights(next_dancer, k) =1/Particles;
              track_state(next_dancer, k) = 1;
        end

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
    end
    if next_dancer > Dancers
        break
    end
end
hold off

pause(1)

for t = 2 : Frames
    time_step = 0;
    Imwork = imread(strcat('DATA1/frame', int2str(t + 109), '.jpg'),'jpg');
    
    imshow(Imwork)
    
    hold on
    
    candidate_regions = getRegion( double(Imwork), 10, 0);
    centres = zeros(size(candidate_regions,1),2);
    
    for r = 1 : size(candidate_regions,1)
        centres(r,:) = [candidate_regions(r).Centroid(1), candidate_regions(r).Centroid(2)];
    end
    
    colours = zeros(size(candidate_regions,1),3);
    
    redChannel = int64(Imwork(:, :, 1));
    greenChannel = int64(Imwork(:, :, 2));
    blueChannel = int64(Imwork(:, :, 3));
    
    for r = 1 : size(candidate_regions,1)
    	col = mask == r; % for example 2 or 3 or whatever region you want.
        colours(r,:) = [mean(redChannel(col)), mean(greenChannel(col)), mean(blueChannel(col))]; 
    end
    
    while ~all(matched)
        time_step = time_step + 1;
        matched = zeros(Dancers);
        for tracker = 1:Dancers
            new_state = track_pos(tracker,1,t-1,:)
            % [new_state, P, new_weights] = condense_function(track_pos());
            distances = pdist2(centres, [new_state(1), new_state(2)]);
            if any(distances < 100)
                col_shift = pdist2(colours, avg_colour(tracker,:));
                col_shift(distances > 100) = Inf;
                [~,choice] = min(col_shift);
                matched(tracker) = choice;
                %TODO update
                track_pos(tracker, 1,t,:) = [candidate_regions(choice).Centroid(1), candidate_regions(choice).Centroid(2),0,0,0,0];
            else
                disp('FAIL')
                plot(new_state(1) + c, new_state(2) + r,'g.')
            end
        end
    end
    for tracker = 1:Dancers 
        radius = sqrt(candidate_regions(matched(tracker)).Area / pi);
%         for k = 1 : Particles
%             plot(track_pos(tracker,k,1,1), track_pos(tracker, k,1,2),'r.')
%         end
        for c = -0.99*radius: radius/10 : 0.99*radius
          r = sqrt(radius^2-c^2);
          plot(candidate_regions(matched(tracker)).Centroid(1) + c, candidate_regions(matched(tracker)).Centroid(2) + r,'r.')
          plot(candidate_regions(matched(tracker)).Centroid(1) + c, candidate_regions(matched(tracker)).Centroid(2) - r,'r.')
        end
    end
    hold off
    pause(1)
end
