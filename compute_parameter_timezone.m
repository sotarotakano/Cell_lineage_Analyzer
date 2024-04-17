function [mdl_dt_s,mdy_dt_s] = compute_parameter_timezone(schnitzcells, thislineage, timecourse, timezone, delta_fr, delta_t,micronsperpixel)
    mdl_dt_s = zeros(1,length(timezone)-1);
    mdy_dt_s = zeros(1,length(timezone)-1);
    
    for k = 2:length(timezone)
        g_Lfit_total =[]; GFP_fit_total = [];
        % Determining the start schnitz for calculating the elongation rate
        for m = 1:length(thislineage)
            interestone = thislineage(m);
            start_fr = min(schnitzcells(interestone).frames);
            end_fr = max(schnitzcells(interestone).frames);
            time = timecourse(start_fr:end_fr);
            start_time = time(1);
            if start_time <= timezone(k-1)
                break;
            end
        end
        
        done = 0;
        % Keep calculating the elongation rate to its ancestors until the birth time of the
        % schnitz reached to the beginning of each time zone
        while ~done
            interestone = thislineage(m);
            start_fr = min(schnitzcells(interestone).frames);
            end_fr = max(schnitzcells(interestone).frames);
            time = timecourse(start_fr:end_fr);
            end_time = time(end);
        
            mintime = max(timezone(k-1),time(1));
            begin_fr_k = max(find(time <= mintime));
            maxtime = min(timezone(k),time(end));
            final_fr_k= max(find(time <= maxtime));
            
            if ~isempty(begin_fr_k) && ~isempty(final_fr_k)
                n_size = schnitzcells(interestone).SZ/micronsperpixel^2;
                TYs = schnitzcells(interestone).SZ.*schnitzcells(interestone).MYs; 
            
                n_size_k = n_size(begin_fr_k:final_fr_k);
                GFP_change_k = TYs(begin_fr_k:final_fr_k);
                time_k = time(begin_fr_k:final_fr_k) - time(begin_fr_k);
        
                [mdl_dt, mdy_dt, g_Lfit_k,GFP_fit_k,~] = compute_elongation_GFP_log_movmean2(n_size_k, GFP_change_k, time_k, delta_fr, delta_t);
                g_Lfit_total = [g_Lfit_total ones(1,length(n_size_k))*mdl_dt];
                GFP_fit_total = [GFP_fit_total ones(1,length(GFP_change_k))*mdy_dt];
            end
            
            if end_time >= timezone(k)
                done = 1;
            end
        
            % Go forward on the lineage to the offspring schnitz
            if m > 1
                m = m - 1;
            else
                break;
            end
        
        end
        
        g_Lfit_total = rmmissing(g_Lfit_total);
        GFP_fit_total = rmmissing(GFP_fit_total);
        if ~isempty(g_Lfit_total)
            mdl_dt_s(k-1) = mean(g_Lfit_total);
            mdy_dt_s(k-1) = mean(GFP_fit_total);
        else
            mdl_dt_s(k-1) = NaN;
            mdy_dt_s(k-1) = NaN;
        end
    
    end
end
