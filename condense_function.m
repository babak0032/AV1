function [new_state, new_P, new_weights, new_x] = condense_function(state, P, weights, observation, dt)
    % Condenstation tracking: Given the previos states and errors, predict the
    % new ones
    % x: State Vector (100, 4)
    % OP: (Old) Covariance of state vec (100,4,4)
    % tracks: Which state are we in (100,1) #TODO Might add later
    % weights: probability of the particles (100,1)
    % z: Observation (3,1) (x, y, radius)
    % 
    % Returns:
    % new_x: State Vector of next time step (100, 4)
    % P: Covariance of state vec of next time step (100,4,4)
    % new_weights: Probability of the new particles

    % Kalman filter static initializations
    R = [0.4, 0.1;
         0.1, 0.4];
    H = [1, 0, 0, 0;
         0, 1, 0, 0];
    Q = eye(4)*1;
    Q(3,3)= 0.01;
    Q(4,4)= 0.01;

    % Tranisiton Matrix
    A = [1, 0, 1, 0;
       0, 1, 0, 1;
       0, 0, 1, 0 ;
       0, 0, 0, 1];  % For all three movements (Stationary, Rotation, and normal Movement)
    B=[0.5, 0;
       0, 0.5;
       1, 0 ;
       0, 1];

    c = [0,0]';
    
    % Number of particles
    NCON = size(weights, 2);
    
    pstationary = 0.0;      % probability of stopping 
    protate_left = 0.0;    % probability of rotating with a partner
    protate_right = 0.0;    % probability of rotating with a partner
    

    % Init answers
    new_state = state;
    new_P = P;
    new_weights = weights;
    
    
    temp_state = zeros(NCON, 4);
    temp_P = zeros(NCON, 4, 4);
    temp_weights = zeros(NCON, 1);
    for t = 1 : dt
        % first create an auxiliary array ident() containing state vector j
        % SAMPLE*p_k times, where p is the estimated probability of j
        SAMPLE = 10;
        ident = zeros(100 * SAMPLE,1);
        idcount=0;

        for i = 1 : NCON    % generate sampling distribution
            num = floor( SAMPLE * 100 * weights(i));  % number of samples to generate
            if num > 0
                ident(idcount+1 : idcount+num) = i * ones(1,num);
                idcount = idcount+num;
            end
        end

        for i = 1 : NCON
            k = ident(ceil(idcount * rand(1))); % select which old sample
            cur_P = squeeze(new_P(k,:,:));
            cur_particle = new_state(k,:)';
            for n = 1: 4
                cur_particle(n) = cur_particle(n) + 5 * sqrt( cur_P(n,n)) * randn(1);
            end
            % Check for tracks(?), which state are we going to be in:
            % this implementation is without tracking 
            % Update new state vector
            r = rand(1);
            if r < pstationary
                cur_particle(3:4) = 0;
            elseif r < (protate_left + pstationary)
                angle = atan2(cur_particle(4), cur_particle(3)) - pi / 2; % I am not sure about this one; z might not be the center at all;
                acc_x = cos(angle)*5; % accelration in x direction
                acc_y = sin(angle)*5; % accelration in x direction
                c = [acc_x,acc_y]';
            elseif r < (protate_left + protate_right + pstationary)
                angle = atan2(cur_particle(4), cur_particle(3)) + pi / 2; % I am not sure about this one; z might not be the center at all;
                acc_x = cos(angle)*5; % accelration in x direction
                acc_y = sin(angle)*5; % accelration in x direction
                c = [acc_x,acc_y]';
            else % normal motion
            end         

            temp_state(i,:) = A * cur_particle + B * c;      % predict next state vector
            
            cur_P = A * cur_P * A' + Q;    % predicted error
            % corrections
            K = cur_P * H' * inv(H * cur_P * H' + R);      % gain
            temp_state(i,:) = temp_state(i,:) + (K * (observation' - H * temp_state(i,:)'))';    % corrected state
            temp_P(i,:,:) = (eye(4) - K * H) * cur_P;                    % corrected error

            % weight hypothesis by distance from observed data
            dvec = observation - [temp_state(i,1), temp_state(i,2)];
            temp_weights(i) = 1/(dvec*dvec');
        end
        totalw = sum(temp_weights);
        new_state = temp_state;
        new_P = temp_P;
        new_weights = temp_weights/totalw;
    end
    new_x = sum( repmat(new_weights,1,4)' .* new_state',2);
end