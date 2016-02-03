function [new_x, temp_state] = predict(state, weights, dt)
    % Condenstation tracking: Given the previos states and errors, predict the
    % new ones
    % state: State Vector (x,y,vx,vy) (100, 4)
    % P: (Old) Covariance of state vec (100,4,4)
    % weights: probability of the particles (100,1)
    % 
    % Returns:
    % new_x: State Vector of next time step (1, 4)

    % Tranisiton Matrix
    A=[1, 0, 1, 0;
       0, 1, 0, 1;
       0, 0, 1, 0 ;
       0, 0, 0, 1];  % For all three movements (Stationary, Rotation, and normal Movement)
    B=[0.5, 0;
       0, 0.5;
       1, 0 ;
       0, 1];
   
   
    c = [0,0]';
    
    % Number of particles
    NCON = size(weights, 1);

    pstationary = 0.05;      % probability of stopping 
    protate_left = 0.15;    % probability of rotating with a partner
    protate_right = 0.15;    % probability of rotating with a partner
    
    temp_state = state;
    for t = 1:dt
        for i = 1 : NCON
            % Update new state vector
            r = rand(1);
            if r < pstationary
                temp_state(i,3:4) = 0;
            elseif r < (protate_left + pstationary)
                angle = atan2(temp_state(i,4), temp_state(i,3)) + pi / 2; % I am not sure about this one; z might not be the center at all;
                acc_x = cos(angle); %accelration in x direction
                acc_y = sin(angle); %accelration in x direction
                c = [acc_x,acc_y]';
            elseif r < (protate_left + protate_right + pstationary)
                angle = atan2(temp_state(i,4), temp_state(i,3)) - pi / 2; % I am not sure about this one; z might not be the center at all;
                acc_x = cos(angle); %accelration in x direction
                acc_y = sin(angle); %accelration in x direction
                c = [acc_x,acc_y]';
            else % normal motion
            end
            temp_state(i,:) = A * temp_state(i,:)' + B * c;      % predict next state vector
        end
    end
    new_x = sum( repmat(weights,4,1) .* temp_state',2);

end