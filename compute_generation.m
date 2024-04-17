%% Compute generation
    function [thislineage, generation] = compute_generation(schnitzcells, schnitz)
        generation = 1; thislineage = [schnitz]; done = 0;
        while ~done
            if (schnitzcells(schnitz).P ~= 0)
                generation = generation+1;
                thislineage = [thislineage schnitzcells(schnitz).P];
            end
            schnitz = schnitzcells(schnitz).P;
            done = (schnitz <=0);
        end
        
    end