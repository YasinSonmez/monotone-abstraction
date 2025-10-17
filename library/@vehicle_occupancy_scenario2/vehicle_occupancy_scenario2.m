classdef vehicle_occupancy_scenario2 < handle
    
    % Description
    % Class that keeps track of which zones in the intersection are being
    % occupied by a particular vehicle.
    
    properties (SetAccess=protected)
        zone;
    end
    
    methods
    
        function ego_occupancy = vehicle_occupancy_scenario2(s)
        
            % Description
            % return the id of the 'intersection zone' that the vehicle is currently
            % occupying, where zones are defined as follows:
            %
            % - zone 1 := {s: s < con.collision_zone_width}
            % - zone 2 := {s: abs(s) <= con.collision_zone_width}
            % - zone 3 := {s: s > con.collision_zone_width}
            %
            % where s is the current position along the ego vehicle's path, and s = 0
            % marks the point where the vehicle's path crosses with the path of the
            % other vehicle, and con.collision_zone_width is the width of the
            % collision zone (as defined in the paper "Making intersections safer
            %                                             with I2V communication.")
            %
            % intuitive names / descriptions of the zones are given in the file
            %
            %
            % Inputs
            % s - the position of the vehicle along its path
            %
            % Outputs
            % zone_id - the id of the zone occupied by the vehicle, based on the
            %           above description
            %
            
            % use modelling parameters defined in constants file
            con = veh_turn_const_v4;
            
            if(s < -con.collision_zone_width)
                ego_occupancy.zone = 1;
            elseif(abs(s) <= con.collision_zone_width)
                ego_occupancy.zone = 2;
            else
                ego_occupancy.zone = 3;
            end
        end
        
        function is_safe = is_ego_safe(ego_occupancy, oncoming_occupancy)
            % we are unsafe if:
            % 1) the ego vehicle cut in front of the oncoming vehicle, or
            % 2) both vehicles are occupying the intersection simultaneously
            is_safe = ~((ego_occupancy.zone > oncoming_occupancy.zone) || ...
                        (ego_occupancy.zone == 2 && ...
                         oncoming_occupancy.zone == 2));
        end
        
        function reached_goal = has_ego_reached_goal(ego_occupancy, ...
                                                      oncoming_occupancy)
            % ego vehicle successfully reached goal set if:
            % 1) the ego vehicle is safe, and
            % 2) it let the oncoming vehicle reach zone 3
            reached_goal = ...
                ego_occupancy.is_ego_safe(oncoming_occupancy) && ...
                (oncoming_occupancy.zone == 3);
        end
        
    end
end
