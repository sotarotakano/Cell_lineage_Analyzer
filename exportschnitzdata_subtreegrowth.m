function data  = exportschnitzdata_subtreegrowth(schnitzcells,p,whichones,filepath,varargin)
% This function computes parameters regarding the cellular growth 
% based on microscopic images and basically from the view of subtree growth. 
% This program is compatible to the data structure
% of Schnitzcellls (Young et al., 2012 Nat. Protoc.) and 
% should be run in MATLAB environemnt. 
% Users need to prepare a "schnitzcells" struct variable before running
% this function by "compileschnitz" function.

% This function returns the parameters regarding subtree growth in Takano
% et al. 2024. All the initially set parameters are same as those indicated
% in the original paper.

% In addition to the "schnitzcells" struct variable, users need to prepare
% timecourse.csv file (please see a test data).

%% Initial setting of parameters.
% micronsperpixel: Initially, Schnitzcells program uses "pixel" as a unit of length or size of cell. 
%                          Therefore the users should give a coefficient for convert pixel to micrometer.

micronsperpixel = 11.0132;

% delta_fr (frames): The number of frames used for calculating the size increase rate (i.e., used frames for
%                           computing moving mean). In the original paper, differenet
%                           parameters were used between high-nutient (delta_fr_nstv)
%                           and low-nutrient (delta_fr_stv).

delta_fr_nstv = 3;
delta_fr_stv = 9;

% delta_t (mins):  The time interval for computing a size increment (deltaS) of a cell. 
%                        Basically, deltaS is given as  S(t + delta_t) - S(t), and
%                        different delta_t parameters were given in tow different
%                        conditions same as the case of delta_fr.

delta_t_nstv = 30;
delta_t_stv = 400;

% switchtime_0 (mins): Start time of experiments.
% switchtime_1 (mins): Timing from high-nutrient (initial) to low-nutrient.
% switchtime_2 (mins): Timing from low-nutrient to high-nutrient (recovery).
% recovery_period_duration (mins): The duration of the recovery period.
switchtime_0 = 0;
switchtime_1 = 360;
switchtime_2 = 4680;
recovery_period_duration = 600;

analysis_starttime_nstv = 120;
analysis_endtime_stv = 3520;


% timezone (mins): Here, we divided the time into 7 periods as defined by this
%                           parameter. In Takano et al., we first exposed cells in
%                           high-nutrient period for 360 mins, and thus the
%                           first time zone is set as 0 -360. Then,
%                           low-nutrient period started and continued for
%                           4320 mins, and this period is further divided
%                           by 6 periods as shown in an example.

timezone = [360 1080 1800 2520 3240 3960 4680];

% parameter: this is the target for the subtree growth analysis (default).
% Users also can use "GFP", that is the total GFP intensity increase rate
% in each subtree population.
parameter = "growth";
%% Mainbody
% Loading timecourse data from the directory.
dirpath = [filepath '\' p.movieDate '\' p.movieName];
disp(['The information in ' dirpath ' is used for the analysis.']);
timecoursedata = csvread([dirpath '\timecourse.csv']);

% Checking which schnitz(s) are used for the analysis. Here, only approved
% schnitz are subject to the further analysis.
if isempty(whichones)
    whichones = find([schnitzcells.approved]);
end

nonapproved =  find([schnitzcells.approved]==0);
if ~isempty(nonapproved)
    for i = 1:length(nonapproved)
        disp(['Schnitz ' nonapproved(i) ' is not approved yet.']);
    end
end
whichones = sort(whichones);

% Loading timecourse and timeframe data from a given csv file. 
timecourse = timecoursedata(:,1).';
timeframes = timecoursedata(:,2).';
timecourse = [0 timecourse];

% Searching for switch times from timecouse information.
switchtime_1 = max(timecourse(timecourse <= switchtime_1));
switchtime_2 = max(timecourse(timecourse <= switchtime_2));
endtime = switchtime_2 + recovery_period_duration;

% experimental number (colony number)
exp_No = str2double(p.movieName);

% time frame numbers for analyzing elongation rate
dS_dt_number = cell(1,length(timezone)-1);
for i = 1:length(timezone)-1
    dS_dt_number(1,i) = {['dS_dt___' num2str(i)]};
end

VariableNames = ['Exp_No' 'Origin' 'Schnitz' 'Parent' 'Grandparent' 'Sister' 'birth_time' 'duration' ...
    'doubling_time' 'division_time' 'daughter1' 'daughter2' 'lag_time' ...
    'dS_dt_highN' 'dS_dt_lowN' 'dS_dt_recovery' dS_dt_number ...
    'f_offsprings' 't_offsprings' 'offsprings' 'v_offsprings' 'lysed_offsprings'];

output = zeros(length(whichones),length(VariableNames));

for i = 1:length(whichones)
    Origin = NaN; bt = NaN; d = NaN;  daughter1 = NaN; daughter2 = NaN;
    dt = NaN;  grandparent =NaN; parent = NaN; sister = NaN; 
    division_t = NaN; lag_t =NaN; dS_dt_nstv = NaN; dS_dt_stv = NaN; dS_dt_recovery = NaN; 
    mdS_dt_stv = NaN(1,length(timezone)-1); mdeltaL = NaN(1,length(timezone)-1);
    t_offsprings = 0; v_offsprings = 0; offsprings = 0; lysed_offsprings = 0; f_offsprings = 0;
   
    % Most of the parameters, especially those which are independent of
    % lineage and genealogical information are computed by
    % a "schnitzanalysis" function. See below.
    
    schnitzanalysis(whichones(i));
    
    % Adding genealogical information for the target schnitz
    if schnitzcells(whichones(i)).P ~= 0
        parent = schnitzcells(whichones(i)).P;
        grandparent = schnitzcells(parent).P;
        sister = schnitzcells(parent).D;
        if sister == whichones(i)
            sister = schnitzcells(parent).E;
        end
    end
    
    % Computing subtree growth rate: For the schnitz born before the second
    % switchtime (i.e., existed at the end of the low-nutrient period), we
    % computed "subtree size increase rate". Subtree is a set of cell-
    % lineages derived from a common ancestral cell. Here, we computed the
    % increase rate of a total cell area occupied by whole cell population belonging to each subtree. 
    
    % We computed the subtree size increase rate for the cells whose birth
    % time (bt) was earlier than the switchtime_2.
    if bt < switchtime_2
        
        % Computing subtree size increase rate in high-nutrient period.
        % Thus, the start time point for the analysis is either of
        % user defined analysis start time (analysis_starttime_nstv, default: 120 min) 
        % or the birth time of the cell (ancestral cell of a subtree) of interest.
        
        if  bt +d > switchtime_0 &&  bt < switchtime_1
            starttime = max(bt,analysis_starttime_nstv);
            
            % t_len: timecourse information used for the subtree growth analysis 
            t_len = timecourse(max(find(timecourse <= starttime)):max(find(timecourse <= switchtime_1)));
            length_increment = [];
            
            % length_increment: an array storing the timecourse data of areas
            % occupied by a subtree of interest
            for k = 1: length(t_len)
                total_increment = 0;
                total_amount(whichones(i),t_len(k),parameter);
                 if total_increment == 0
                    disp(['The schnitz ' num2str(whichones(i)) ' disappeared (lysed) at t = ' num2str(t_len(k)) ' during its growth...']);
                    t_len = t_len(1:length(length_increment));
                    break;
                end
                length_increment = [length_increment total_increment];
            end
            [dS_dt_nstv, ~, dlogS_dt, ~] = compute_elongation_GFP_log_movmean2(length_increment, length_increment, t_len,delta_fr_nstv,delta_t_nstv);
        end
        
        % Computing subtree size increase rate in low-nutrient period.
        % Thus, the start time point for the analysis is either of
        % switchtime_1 or the birth time of the cell (ancestral cell of a subtree) of interest.
        starttime = max(switchtime_1, bt);
        
        % Here, the used timeframe for the analysis of subtree growth is
        % limited by analysis_endtime_stv (this is set to 3520 min in the
        % original paper), and the subtree growth data from
        % "switchtime_1" to "analysis_endtime_stv" is used.
        
        analysis_endfr = find(timecourse <= analysis_endtime_stv,1,'last') - max(find(timecourse <= starttime));
        
        % t_len: timecourse information used for the subtree growth analysis 
        t_len = timecourse(max(find(timecourse <= starttime)):find(timecourse <= switchtime_2, 1, 'last' ));
        
        % length_increment: an array storing the timecourse data of areas
        % occupied by a subtree of interest
        length_increment = [];
        for k = 1: length(t_len)
            % Here, "total_increment" is defined as the total area occupied
            % by subtree population of interest at time t_len(k)
            total_increment = 0;
            total_amount(whichones(i),t_len(k),parameter);
            if total_increment == 0
                % If the total_increment is estimated as 0, which means
                % that the subtree of interest totally disappeared at that
                % time point (e.g., all the cells were lysed etc ...).
                % Then, the estimation is stopped at that time point and
                % move to the next step.
                disp(['The schnitz ' num2str(whichones(i)) ' disappeared (lysed) during its growth...']);
                t_len = t_len(1:length(length_increment));
                break;
            end
            length_increment = [length_increment total_increment];
        end
        
        % Computing subtree size increase rate in each time zone.
        for k = 1:length(timezone)-1
            start_fr_k = max(find(t_len <= timezone(k)));
            end_fr_k = max(find(t_len <= timezone(k+1)));
            
            if isempty(start_fr_k) || isempty(end_fr_k)
                continue;
            end
            length_increment_k = length_increment(start_fr_k:end_fr_k) ;
            mdeltaL(k) =length_increment_k(end) - length_increment_k(1);
            t_len_k = t_len(start_fr_k:end_fr_k);

            [mdS_dt_stv(k), ~, ~, ~] = compute_elongation_GFP_log_movmean2(length_increment_k, length_increment_k, t_len_k,delta_fr_stv,delta_t_stv);
        end
       
       % The subtree growth data in low-nutrient period is finally trimmed to a given timepoint.
       if analysis_endfr > 0 
            if length(length_increment) >= analysis_endfr
                length_increment_stv = length_increment(1:analysis_endfr);
                t_len_stv = t_len(1:analysis_endfr);
            end
                
            if ~isempty(t_len)
                [dS_dt_stv, ~, dlogS_dt, ~] = compute_elongation_GFP_log_movmean2(length_increment_stv, length_increment_stv, t_len_stv,delta_fr_stv,delta_t_stv);
            end
       end
        
        %  Computing subtree size increase rate in recovery period.
        t_len = timecourse(max(find(timecourse <= switchtime_2)):max(find(timecourse <= endtime)));
        
        length_increment = [];
        for k = 1: length(t_len)
            total_increment = 0;
            total_amount(whichones(i),t_len(k),parameter);
            if total_increment == 0
                disp(['The schnitz ' num2str(whichones(i)) ' disappeared (lysed) during its growth...']);
                t_len = t_len(1:length(length_increment));
                break;
            end
            length_increment = [length_increment total_increment];
        end
        if ~isempty(length_increment)
            [dS_dt_recovery, ~, ~, ~] = compute_elongation_GFP_log_movmean2(length_increment, length_increment, t_len,delta_fr_nstv,delta_t_nstv);
        end   
    end
 
    thisone_output = [exp_No Origin whichones(i) parent grandparent sister bt d ...
        dt division_t daughter1 daughter2 lag_t ...
        dS_dt_nstv dS_dt_stv dS_dt_recovery mdS_dt_stv ...
        f_offsprings t_offsprings offsprings v_offsprings lysed_offsprings];
    
    output(i,:) = thisone_output;
    
end

data = array2table(output);
data.Properties.VariableNames = VariableNames;

for i = 1:height(data)
    t_offsprings = 0; v_offsprings = 0;
    offsprings = 0; lysed_offsprings = 0;
    offsprings_lag = 0; 

    countoffsprings_before_recovery(i);
    countoffsprings_total(i);

    % f_offsprings: the number of final offsprings at the end of experiments
    % offsprings: the number of offsprings at the start of the recovery period
    % t_offsprings: the number of total offsprings produced in the recovery period
    % v_offsprings: the number of viable offsprings (i.e. can grow in the recovery period)
    % lysed_offsprings: the number of offsprings that disappeared during growth
    % lagtime_average: the average of lag time in starter cells belonging to each subtree 
    
    data.lagtime_average(i) = offsprings_lag/offsprings;
    data.f_offsprings(i) = t_offsprings;
    data.t_offsprings(i) = t_offsprings - offsprings;
    data.offsprings(i) = offsprings;
    data.v_offsprings(i) =v_offsprings;
    data.lysed_offsprings(i) =lysed_offsprings;
end

data(data.lagtime_average == 0,'lagtime_average') = {NaN};

    %% function schnitzanalysis processed tracked data to extract the information about cell growth and fluorescence intensity.
    function schnitzanalysis(thisone)
        % the number of frames (length of movie) for this schnitz
        frames = length(schnitzcells(thisone).frames);
        
        % time course of schnitz
        start_fr = min(schnitzcells(thisone).frames);
        end_fr = max(schnitzcells(thisone).frames);
        time = timecourse(start_fr:end_fr);
        
        % cell length
        n_size = schnitzcells(thisone).SZ/micronsperpixel^2;
        
        % duration
        d = timecourse(end_fr) - timecourse(start_fr);
        
        % Birth time
        bt = timecourse(start_fr);
        
        % YFP and RFP signal changes (total fluorescence per cell)
        TYs = schnitzcells(thisone).SZ.*schnitzcells(thisone).MYs;
        
        % Distinguishing whether each schnitz is dautgherless or no
        if (schnitzcells(thisone).D ~= 0)
            D = 1;
            daughter1 = schnitzcells(thisone).D;
            daughter2 = schnitzcells(thisone).E;
            d = d + timeframes(end_fr);
            %division time: the time of division event
            division_t = bt + d;
            dt = d;
        else
            D = 0;
            division_t = NaN;
            dt = NaN;
            daughter1 = NaN;
            daughter2 = NaN;
        end
        
        if bt == 0
            dt = NaN;
        end
        
        % parent and grandparent
         if schnitzcells(thisone).P ~= 0
            parent = schnitzcells(thisone).P;
            grandparent = schnitzcells(parent).P;
         end
        
        % lagtime: the duration from the nutrient addition to first division
        % This value is derived for the schnitz existing after the swiching time
        if bt < switchtime_2 && bt +d >= switchtime_2
            lag_t = min(endtime,division_t)-switchtime_2;
        else
            lag_t = NaN;
        end
    end

    %% Counting the number of offsprings before entrering recovery period
    function countoffsprings_before_recovery(daughterone)
        last_time = data.birth_time(data.Schnitz == daughterone) + data.duration(data.Schnitz == daughterone);
        if last_time<=  switchtime_2
            if  ~isnan(data.daughter1(data.Schnitz == daughterone))
                countoffsprings_before_recovery(data.daughter1(data.Schnitz == daughterone));
                countoffsprings_before_recovery(data.daughter2(data.Schnitz == daughterone));
            end
        else
            offsprings = offsprings+1;
            
            if ~isnan(data.lag_time(data.Schnitz == daughterone))
                offsprings_lag = offsprings_lag + data.lag_time(data.Schnitz == daughterone);
            end
            
           if  data.division_time(data.Schnitz == daughterone) < endtime && data.birth_time(data.Schnitz == daughterone) < switchtime_2
                v_offsprings = v_offsprings + 1;
           end
        end
    end
    %% Counting the number of offsprings (total)
    function countoffsprings_total(daughterone)
        last_time = data.birth_time(data.Schnitz == daughterone) + data.duration(data.Schnitz == daughterone);
        if  ~isnan(data.daughter1(data.Schnitz == daughterone)) && data.division_time(data.Schnitz == daughterone) <= endtime
            countoffsprings_total(data.daughter1(data.Schnitz == daughterone));
            countoffsprings_total(data.daughter2(data.Schnitz == daughterone));
        else
           if last_time < endtime && data.birth_time(data.Schnitz == daughterone) < switchtime_2
                lysed_offsprings  = lysed_offsprings + 1;
            else
                t_offsprings = t_offsprings+1;
           end
        end
    end
    %% total cell length at time t
    function total_amount(schnitznum,timepoint,params)
        start_fr = min(schnitzcells(schnitznum).frames);
        end_fr = max(schnitzcells(schnitznum).frames);
        
        if params == "growth"
            target = schnitzcells(schnitznum).SZ/micronsperpixel^2;
        elseif params == "GFP"
            target = schnitzcells(schnitznum).SZ.*schnitzcells(schnitznum).MYs;
        end
        timepoint = max(timecourse(timecourse <= timepoint));
        
        time = timecourse(start_fr:end_fr);
        
        if time(end) >= timepoint && time(1) <= timepoint
            total_increment = total_increment +  target(find(time == timepoint));
            
        elseif schnitzcells(schnitznum).D ~= 0
            total_amount(schnitzcells(schnitznum).D,timepoint,params);
            total_amount(schnitzcells(schnitznum).E,timepoint,params);
        end
    end

end