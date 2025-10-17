function v = outer_approx_boundary(h, v_L, a_min_brake, a_max_brake, d_min)
    
% Consider a vehicle-following scenario with an ego vehicle and a lead
% vehicle directly ahead of it. We can model the dynamics as:

% \dot{h}   = v_L - v,
% \dot{v}   = a,
% \dot{v_L} = a_L,

% where h = headway, v = ego velocity, v_L = lead velocity,
% a = ego acceleration, and a_L = lead acceleration.

% Given values for h and v_L, this function returns the value for v such 
% that the point (h, v, v_L) lies on the boundary of an over-approximation
% of the set of safe states in this scenario.

% See "Improving Urban Traffic Throughput With Vehicle Platooning:
%      Theory and Experiments" for the definition of the safe set.

% Check inputs
if(h < d_min)
    error(['Input h = ', h, ' is below the minimum distance d_min = ', d_min, '.']);
end
if(v_L < 0)
    error(['Input v_L = ', v_L, ' is negative. ' ...
           'Only nonnegative values are permitted.']);
end
if(a_min_brake < 0)
    error(['Input a_min_brake = ', a_min_brake, ' is negative. ' ...
           'Only nonnegative values are permitted.']);
end
if(a_max_brake < 0)
    error(['Input a_max_brake = ', a_max_brake, ' is negative. ' ...
           'Only nonnegative values are permitted.']);
end

% Return result based on safe set over-approximation definition.
v_squared = 2*a_max_brake*(h + v_L^2/(2*a_min_brake) - d_min);
if(v_squared < 0)
    v = 0; % this state is unsafe
else
    v = sqrt(v_squared);
end

end
