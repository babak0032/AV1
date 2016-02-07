function answer = evaluation(observations, path)
% Plot every thing
% observations: all 4 observation (They should be in the correct order, but
% it doesn't matter) (region, frame, [x,y]) : (4, 210, 2)
% path:(dancers, frames, [x,y]) : (4,210,2)

% vector of number of correct detetctions for each dancer

% dist_stats = zeros(4,size(observations,2));

% Number of detections withing every ground truth pixels:
num_of_detections = zeros(1,4);

% Every time we failed to detect anything within 10 pixels of ground truth
% pixels
incorrect_matchings = zeros(1,4);
incorrect_tracking = zeros(1,4);

load('DATA1/positions1.mat')
answer = 1;
index = 1;

% plot the ground truth positions (4,210,2)
% plot(squeeze(positions(1,1:209,1)), squeeze(positions(1,1:209,2)))
% hold on
% plot(squeeze(observations(1,:,1)), squeeze(observations(1,:,2)))
% pause(2)

for frame=1:size(observations,2)  
  for dancer = 1:4
      actual_pos = squeeze(positions(dancer,frame,:));
      min_dist = inf;
      for detection =1:4
          dist = pdist2(actual_pos', squeeze(observations(detection,frame,:))', 'euclidean');
          
          if dist < 10
              num_of_detections(dancer) = num_of_detections(dancer) + 1;
              dist_data(index) = dist;
              index = index + 1;
          end
          
          if dist < min_dist
              min_dist = dist;
          end          
      end
      if min_dist > 10
         incorrect_matchings(dancer) = incorrect_matchings(dancer) + 1;
         min_dist;
         %disp(frame+109)
      end
      
      % dist_stats(dancer,frame) = min_dist;      
  end  
  
  actual_pos = squeeze(positions(:,frame,:));
  
  % the dancer are not in the right order, sowe have to re arange them
  dancers_pos = squeeze(path(:,frame,:));
  dancers_pos([2,4],:) = dancers_pos([4,2],:);
  dancers_pos([3,4],:) = dancers_pos([4,3],:);
 

  % Their dancer are in differenet order
  distances = pdist2(actual_pos, dancers_pos, 'euclidean');
  size(distances); % (4,4)
  
  actual_distances = min(distances,[],2)';
  tracker_distances = diag(distances)';
  
  for tracker=1:4

    if actual_distances(tracker) ~= tracker_distances(tracker)
        distances;
        actual_distances;
        tracker_distances;
        incorrect_tracking(tracker) = incorrect_tracking(tracker) + 1;
    end  
  end  
end


num_of_detections
incorrect_matchings
incorrect_tracking

disp('detection error')
sum(incorrect_matchings)/(4*210)
1 - sum(incorrect_matchings)/(4*210)
disp('tracking error')
1 - sum(incorrect_tracking)/(4*210)

size(dist_data)
mean(dist_data)
var(dist_data)
median(dist_data)

end