function [mdS_dt, mdy_dt, delta_S, delta_TY, time_change] = compute_elongation_GFP_log_movmean2(length_change,GFP_change, time_change, delta_fr, delta_t)

%min_frames = delta_fr*2 + 1; 
delta_S = []; delta_TY = []; dS = []; dTY = [];
done = 0;

% First calculated deltaS in a given time frame
if time_change(end)-time_change(1) > delta_t
    % Here, we calculated the size increase rate to (a final frame - delta_fr)
    i = 1;
    while ~done
        t = time_change(i);
        t_1 = t + delta_t;
        i_1 = find(time_change >= t_1, 1);
        
        if isempty(i_1)
            break;
        end
        
        dS = [dS log(length_change(i_1)/length_change(i))/(time_change(i_1)-time_change(i))];
        dTY =[dTY log(GFP_change(i_1)/GFP_change(i))/(time_change(i_1)-time_change(i))];
        
        if i_1 == length(time_change)
            break;
        else
             i = i + 1;
        end
    end
    
    % Then calculating moving mean of each parameters
    delta_S = movmean(dS,delta_fr);
    delta_TY =movmean(dTY,delta_fr);
    
    % if delta_S is negative, it's regarded as 0.
    delta_S(delta_S < 0) = 0;

    mdS_dt = mean(delta_S); mdy_dt = mean(delta_TY);
    
elseif length(length_change) >= 3
    mdS_dt = max(0,log(length_change(end)/length_change(1))/(time_change(end)-time_change(1)));
    mdy_dt = log(GFP_change(end)/GFP_change(1))/(time_change(end)-time_change(1));
else
    mdS_dt = NaN; mdy_dt =NaN;
end

end