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
    if ~any(pdist2(Ignore_pos, repmat([candidate_regions(region).Centroid(1), candidate_regions(region).Centroid(2)], size(Ignore_pos,2),1)) < 50)

        % TODO        avg_colour(next_dancer) = mean( pixels );
        
        avg_size(next_dancer) = candidate_regions(region).Area;
        for k = 1 : Particles
              track_pos(next_dancer, k,1,:) = [floor(candidate_regions(region).Centroid(1) + 15*rand(1)), floor( candidate_regions(region).Centroid(2) + 15*rand(1)),0,0,0,0];
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
% 
% % for t = 2 : Frames
% %     time = 0;
% %     candidate_regions = getRegion( double(imread(strcat('DATA1/frame', int2str(t + 109), '.jpg'),'jpg')), 10, 0);
% %     while ~all(matched)
% %         time = time + 1;
% %         matched = zeros(Dancers);
% %         for tracker = 1:Dancers
% %             if ~matched(tracker)
% %                 prediction = condense_function(track_pos())
% %                 for r = 1 :size(candidate_regions)
% %                     if tests
% % 
% %                         matched(tracker) = 1;
% %                         break
% %                     end
% %                     if matched(tracker)
% %                         updatePARTICLES
% %                     else
% %                         break
% %                     end
% %                 end
% %             end
% %         end
% %     end
% %     drawFindings
% % end
