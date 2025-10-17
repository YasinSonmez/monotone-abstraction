function unsafe = ...
    unsafe_collision_check(h, v, v_L, a_lead_brake, a_ego_brake, v_allow, t)
% Determines if an unsafe collision can occur, given the initial states
% at time t = 0: h = headway, v = follower velocity, v_L = lead velocity.

% Here, a_lead_brake is the deceleration rate of the lead vehicle while 
% braking, and a_ego_brake is the deceleration rate of the ego vehicle
% while braking. Furthermore, v_allow is the
% maximum allowable relative velocity at the time of impact - that is, 
% we require v(t) - v_L(t) <= v_allow for a collision at time t to be
% deemed safe.

% if t = -1, we check if an unsafe collision occured for time t <= 0
% if t = 1, we check if an unsafe collision can occur for time t >= 0
% if t = 0, we ignore vehicle dynamics and check if 
%           h <= 0 and v - v_L <= v_allow

% for the formulas below, we use |a_ego_brake| and |a_lead_brake|
a_ego_brake = abs(a_ego_brake);
a_lead_brake = abs(a_lead_brake);

% for floating point comparisons
threshold = 1e-5;

if(t == -1 || t == 1)
    % we are in a corner case if h = 0
    if(abs(h) < threshold)
        unsafe = (v - v_L > v_allow + threshold);
        return;
    end
    % check inputs
    if(t == 1 && h < 0)
        error('If t = 1, this function assumes a safe initial state.');
    end
    % break into each scenario
    if(t == -1)
        if(h > 0)
            unsafe = false;
            return;
        end
        % check if an unsafe collision can occur for t <= 0
        r = roots([(a_ego_brake - a_lead_brake) / 2, v_L - v, h]);
        rneg = r(r <= 0);
        if(isempty(rneg))
            unsafe = false;
            return;
        end
        t_impact = max(rneg);
    else
        % check if an unsafe collision can occur for t >= 0
        r = roots([(a_ego_brake - a_lead_brake) / 2, v_L - v, h]);
        rpos = r(r >= 0);
        if(isempty(rpos))
            unsafe = false;
            return;
        end
        t_impact = min(rpos);
    end
    % check relative velocity at impact time
    delta_v = ...
            (v - a_ego_brake * t_impact) - (v_L - a_lead_brake * t_impact);
    unsafe = (delta_v > v_allow + threshold);
elseif(t == 0)
    unsafe = ((h <= 0) && (v - v_L <= v_allow));
end

end
