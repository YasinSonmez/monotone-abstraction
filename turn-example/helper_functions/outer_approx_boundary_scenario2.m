function v = outer_approx_boundary_scenario2(s, s0)

% Consider an unprotected left turn scenario where the ego vehicle is 
% attempting to execute a turn at an intersection with an oncoming
% vehicle approaching from the oppposite side of the intersection. We
% can model the system dynamics as:

% \dot{s} = v,
% \dot{v} = a,
% \dot{s0} = v0,

% where s = ego vehicle position along its pathway, v = ego vehicle
% velocity along its pathway, s0 = oncoming vehicle position along its
% pathway.

% Suppose we want the ego vehicle to complete its turn
% (s > con.collision_zone_width) after the oncoming vehicle exits the
% intersection (abs(s0) <= con.collision_zone_width).

% compute safe starting velocity for the ego vehicle, assuming it 
% accelerates at the maximum rate for t >= 0
con = veh_turn_const_v4();
t0 = (con.collision_zone_width - s0) / con.v0_min;

if(t0 < 0)
    % oncoming vehicle already exited the intersection => set v = v_max
    v = con.v_max;
else
    v = (-con.collision_zone_width - s - 0.5*con.a_min_ego*(t0)^2) / t0;
    % do thresholding in case v is out of bounds (not sure if needed)
    v = min(max(v, con.v_min), con.v_max);
end

end
