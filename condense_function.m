function [new_x, P, new_weights] = condense_function(x,OP,weights,z)
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
R=[[0.2845,0.0045]',[0.0045,0.0455]'];
H=[[1,0]',[0,1]',[0,0]',[0,0]'];
Q=0.01*eye(4);
dt=1;

% Tranisiton Matrix
A=[[1,0,0,0]',[0,1,0,0]',[dt,0,1,0]',[0,dt,0,1]'];  % For all three movements (Stationary, Rotation, and normal Movement)
Bu = [0,0,0,0]'; %TODO: It gets updated later with the velocity

% Number of particles/Hypothesis
NCON=100;

pstationary=0.05;      % probability of stopping 
protate=0.30;    % probability of rotating with a partner
xc=zeros(4,1);   % selected state
TP=zeros(4,4);   % predicted covariance


% Init answers
new_x = zeros(NCON,4);
P = zeros(NCON,4,4);
new_weights = zeros(NCON,1);

% first create an auxiliary array ident() containing state vector j
% SAMPLE*p_k times, where p is the estimated probability of j
SAMPLE=10;
ident=zeros(100*SAMPLE,1);
idcount=0;
for i = 1 : NCON    % generate sampling distribution
  num=floor(SAMPLE*100*weights(i));  % number of samples to generate
  if num > 0
    ident(idcount+1:idcount+num) = i*ones(1,num);
    idcount=idcount+num;
  end
end

for i = 1 : NCON
    
    k = ident(ceil(idcount*rand(1))); % select which old sample
    xc = (x(k,:))';
    
    % Update new state vector
    for n = 1 : 4
      xc(n) = xc(n) + 5*sqrt(OP(i,n,n))*randn(1); %TODO: Might need to change this abit
    end
    
    % Check for tracks(?), which state are we going to be in:
    % this implementation is without tracking 
    r = rand(1);
    if r < pstationary
      xc(3) = 0;
      xc(4) = 0;
    elseif r < (protate + pstationary)
      angle = atan2(xc(2)-z(2), xc(1)-z(1)); % I am not sure about this one; z might not be the center at all;
      rad = sqrt((z(2)-xc(2))^2 + (z(1) - xc(1))^2); % (radius) distance from the state/particle to the center of rotation 
      abs_velocity = xc(3)^2 + xc(4)^2;
      acc_x = cos(angle)*(abs_velocity/rad); %accelration in x direction
      acc_y = sin(angle)*(abs_velocity/rad); %accelration in x direction
      Bu = [0,0,acc_x,acc_y]';

    else % normal motion
      continue
    end
    
    % size(A)
    % size(xc)
    % size(Bu)
    
    xp=A*xc + Bu;      % predict next state vector
    TP = reshape(OP(i,:,:), [4,4]);
    
    % size(Q)
    % size(A)
    % size(TP)
    
    PP = A*TP*A' + Q;    % predicted error
    % corrections
    K = PP*H'*inv(H*PP*H'+R);      % gain
    x(i,:) = (xp + K*([z(1),z(2)]' - H*xp))';    % corrected state
    P(i,:,:) = (eye(4)-K*H)*PP;                    % corrected error

    % weight hypothesis by distance from observed data
    dvec = [z(1),z(2)] - [x(i,1),x(i,2)];
    new_weights(i) = 1/(dvec*dvec');
    
    % rescale new hypothesis weights
    totalw=sum(new_weights);
    new_weights=new_weights/totalw;
        
        
end
    
end