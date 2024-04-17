function data  = exportschnitzdata_lineagegrowth(schnitzcells,p,whichones,filepath,varargin)
% function data  = exportschnitzdata_lineagegrowth(schnitzcells,p,whichones,filepath,varargin)
% modified 24/04/16 
% This function computes parameters regarding the cellular growth 
% based on microscopic images. This program is specifically for the data
% analyzed by Schnitzcellls (Young et al., 2012 Nat. Protoc.) and 
% can work in MATLAB environemnt. 
% Users need to prepare a "schnitzcells" struct variable before running
% this function by "compileschnitz" function, which is implemented in the
% original Schnitzcells function. 

% This function returns the most of the growth-related parameters in Takano
% et al. 2024. All the initially set parameters are same as those indicated
% in the original paper.

% In addition to the "schnitzcells" struct variable, users need to prepare
% a timecourse.csv file (please see a test data).

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

% switchtime1 (mins): Timing from high-nutrient (initial) to low-nutrient.
% switchtime2 (mins): Timing from low-nutrient to high-nutrient (recovery).
% recovery_period_duration (mins): The duration of the recovery period.
switchtime_1 = 360;
switchtime_2 = 4680;
recovery_period_duration = 600;

% timezone (mins): Here, we divided the time into 7 periods as defined by this
%                           parameter. In Takano et al., we first exposed cells in
%                           high-nutrient period for 360 mins, and thus the
%                           first time zone is set as 0 -360. Then,
%                           low-nutrient period started and continued for
%                           4320 mins, and this period is further divided
%                           by 6 periods as shown in an example.

timezone = [0 360 1080 1800 2520 3240 3960 4680];

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

% time frame numbers for analyzing the size increae rate
size_increase_rate_number = cell(1,length(timezone)-1);
for i = 1:length(timezone)-1
    size_increase_rate_number(1,i) = {['SizeIncreaseRate__' num2str(i)]};
end

VariableNames = ['Exp_No' 'Origin' 'Schnitz' 'Parent' 'Grandparent' 'Sister' 'birth_time' 'duration' ... 
    'doubling_time' 'division_time' 'daughter1' 'daughter2' 'cell_size_stv_t0' 'generation_high' ...
    'generation_low' 'generation_recovery' 'len_start' 'len_end' 'lag_time' 'size_increase_rate_high' ... 
    'size_increase_rate_low' 'dGFP_dt_high' 'dGFP_dt_low' size_increase_rate_number 'f_offsprings' ...
    't_offsprings' 'offsprings' 'v_offsprings' 'lysed_offsprings' 'starvation_entry_time' 'first_divisiontime_stv'];

output = zeros(length(whichones),length(VariableNames));

for i = 1:length(whichones)
    % initialize parameters 
    Origin = NaN; bt = NaN; d = NaN;  daughter1 = NaN; daughter2 = NaN;
    dt = NaN;  size_start_stv = NaN; grandparent =NaN; parent = NaN; sister = NaN; 
    Gg = NaN; Gr = NaN; Gs = NaN;  l_start = NaN; l = NaN; division_t = NaN; lag_t =NaN;
    mdl_dt_s = NaN(1,length(timezone)-1); dl_dt_g =NaN; dl_dt_s = NaN; dy_dt_g = NaN; dy_dt_s = NaN;
    t_offsprings = 0; v_offsprings = 0; offsprings = 0; lysed_offsprings = 0;  
    starvation_entry_time = NaN;  first_divisiontime_stv = NaN; f_offsprings = 0;
   
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
    
    thisone_output = [exp_No Origin whichones(i) parent grandparent sister bt d dt division_t ...
        daughter1 daughter2  size_start_stv Gg Gs Gr l_start l lag_t dl_dt_g dl_dt_s dy_dt_g dy_dt_s ... 
        mdl_dt_s f_offsprings t_offsprings offsprings v_offsprings lysed_offsprings starvation_entry_time first_divisiontime_stv];
    
    output(i,:) = thisone_output;
    
end

data = array2table(output);
data.Properties.VariableNames = VariableNames;

%% Computing total offsprings in each schnitz
% Computing total offspring cells produced in high- (initial), 
% low-nutrient, or recovery period.
 
for i = 1:height(data)
    t_offsprings = 0; v_offsprings = 0;
    offsprings = 0; lysed_offsprings = 0;
 
    countoffsprings_before_recovery(i);
    countoffsprings_total(i);
    
    % f_offsprings: the number of final offsprings at the end of experiments
    % offsprings: the number of offsprings at the start of the recovery period
    % t_offsprings: the number of total offsprings produced in the recovery period
    % v_offsprings: the number of viable offsprings (i.e. can grow in the recovery period)
    % lysed_offsprings; the number of offsprings that disappeared during growth
    
    data.f_offsprings(i) = t_offsprings;
    data.t_offsprings(i) = t_offsprings - offsprings;
    data.offsprings(i) = offsprings;
    data.v_offsprings(i) =v_offsprings;
    data.lysed_offsprings(i) =lysed_offsprings;
end


    %% function schnitzanalysis processed tracked data to extract the information about cell growth and fluorescence intensity.
    function schnitzanalysis(thisone)
        % the number of frames (length of movie) for this schnitz
        frames = length(schnitzcells(thisone).frames);
        
        % time course of schnitz
        start_fr = min(schnitzcells(thisone).frames);
        end_fr = max(schnitzcells(thisone).frames);
        time = timecourse(start_fr:end_fr);
        
        % cell size
        n_size = schnitzcells(thisone).SZ/micronsperpixel^2;
        l_start = n_size(1);
        l = n_size(end);
        %dl = n_size(end) - n_size(1);
        
        % duration from the birth to the end 
        d = timecourse(end_fr) - timecourse(start_fr);
        
        % Birth time
        bt = timecourse(start_fr);
        
        % YFP or RFP signal changes (total fluorescence per cell)
        TYs = schnitzcells(thisone).SZ.*schnitzcells(thisone).MYs;
        
        %Optional
        %TRs = schnitzcells(thisone).SZ.*schnitzcells(thisone).MRs;
        %MYs = schnitzcells(thisone).MYs;
        %MRs = schnitzcells(thisone).MRs;
        %y = max(schnitzcells(thisone).MYs);                %Mean YFP signals
        %r = max(schnitzcells(thisone).MRs);     %Mean RFP signal
        
        
        % Distinguishing whether each schnitz is dautgherless or no
        if (schnitzcells(thisone).D ~= 0)
            %D = 1;
            daughter1 = schnitzcells(thisone).D;
            daughter2 = schnitzcells(thisone).E;
            d = d + timeframes(end_fr);
            %division time: the time of division event
            division_t = bt + d;
            dt = d;
        else
            %D = 0;
            division_t = NaN;
            dt = NaN;
            daughter1 = NaN;
            daughter2 = NaN;
        end
        
        % We don't consider the case for the schnitz which existed before
        % the start of experiments.
        if bt == 0
            dt = NaN;
        end
        
        % A frame where the transition from high- to low-nutrient occured
        switch_fr_1 = find(time == switchtime_1);
        
         % Finding a parent and a grandparent schnitz
         if schnitzcells(thisone).P ~= 0
            parent = schnitzcells(thisone).P;
            grandparent = schnitzcells(parent).P;
         end
         
        % A frame where the transition to recovery period occured
        switch_fr_2 = find(time == switchtime_2);
        
        % lag_t: a lag time, the duration from the nutrient addition to first
        % division in the recovery period.
        if bt < switchtime_2 && bt +d >= switchtime_2
            lag_t = min(endtime,division_t)-switchtime_2;
        else
            lag_t = NaN;
        end
        
        % Compute generation and the most ancestral schnitz (Origin)
        [thislineage, ~] = compute_generation(schnitzcells,thisone);
        Origin= thislineage(end);
        
        % Computing size increase rate and GFP production rate in
        % high-nutrient condition (initial growth)
        
        % dl_dt_g, dy_dt_g: the size increase rate or GFP production rate,
        % respectively, we computed those paramters for each schnitz, from
        % their birth to division (or the end of the period).
        if  bt < switchtime_1 && frames > 1
            if ~isempty(switch_fr_1)
                [dl_dt_g, dy_dt_g, ~, ~] = compute_elongation_GFP_log_movmean2(n_size(1:switch_fr_1), ...
                    TYs(1:switch_fr_1), time(1:switch_fr_1)-time(1),delta_fr_nstv,delta_t_nstv);
            else
                [dl_dt_g, dy_dt_g, ~, ~] = compute_elongation_GFP_log_movmean2(n_size, TYs, time-time(1), ...
                    delta_fr_nstv,delta_t_nstv);
            end
        else
            dl_dt_g = NaN; dy_dt_g = NaN; 
        end
        
        % Computing elongation rate and GFP production rate
        % after the transition to the low-nutrient condition
        
        % dl_dt_s, dy_dt_s: the size increase rate or GFP production rate,
        % respectively, we computed those paramters for each schnitz, from
        % their birth to division (or the end of the period).
        
        g_Lfit = [];  GFP_fit = [];
        if bt+d > switchtime_1 && bt < switchtime_2
            if  ~isempty(switch_fr_1)
                if  ~isempty(switch_fr_2)
                    [dl_dt_s, dy_dt_s, g_Lfit, GFP_fit, timechange] = compute_elongation_GFP_log_movmean2(n_size(switch_fr_1:switch_fr_2), TYs(switch_fr_1:switch_fr_2), ...
                        time(switch_fr_1:switch_fr_2)-time(switch_fr_1),delta_fr_stv,delta_t_stv);
                else
                    [dl_dt_s, dy_dt_s, g_Lfit, GFP_fit,timechange] = compute_elongation_GFP_log_movmean2(n_size(switch_fr_1:end),TYs(switch_fr_1:end), ...
                        time(switch_fr_1:end)-time(switch_fr_1),delta_fr_stv,delta_t_stv);
                end
            else
                if  ~isempty(switch_fr_2)
                    [dl_dt_s, dy_dt_s, g_Lfit, GFP_fit,timechange] = compute_elongation_GFP_log_movmean2(n_size(1:switch_fr_2), TYs(1:switch_fr_2), ...
                        time(1:switch_fr_2)-time(1),delta_fr_stv,delta_t_stv);
                else
                    [dl_dt_s, dy_dt_s, g_Lfit, GFP_fit,timechange] = compute_elongation_GFP_log_movmean2(n_size,TYs,time-time(1),delta_fr_stv,delta_t_stv);
                end
            end
        else
            dl_dt_s = NaN; dy_dt_s = NaN;
        end
        
        % "Retrospective analysis"
        % From here, we computed the paramters as a cell lineage by tracing
        % back the lineage history of each schnitz.
        
        % Compute generations (i.e., the number of divisions which the cell of interest
        % encounter) during a target condition.
        
        % Gg: generations counted from the start of the high-nutrient period
        Gg = 0;
        for  j = 1:length(thislineage)
            bt_start = timecourse(min(schnitzcells(thislineage(j)).frames));
            if (bt_start <=  switchtime_1) && bt_start > 0
                Gg = Gg + 1;
            end
        end
        
        % Compute generation and length right after the start of starvation
        % Gs: generations counted from the start of the low-nutrient period
        if (bt + d > switchtime_1)
            % initialize Gs
            Gs = 0;
            
            for  j = 1:length(thislineage)
                bt_start = timecourse(min(schnitzcells(thislineage(j)).frames));
                if (bt_start <=  switchtime_1)
                    st_fr = schnitzcells(thislineage(j)).frames(1);
                    en_fr = schnitzcells(thislineage(j)).frames(end);
                    sw_fr = find(timecourse(st_fr:en_fr) == switchtime_1);
                    if ~isempty(sw_fr)
                        size_start_stv = schnitzcells(thislineage(j)).SZ(sw_fr)/micronsperpixel^2;
                        starvation_entry_time = switchtime_1 - timecourse(st_fr);
                        first_divisiontime_stv = timecourse(en_fr) - switchtime_1;
                    else
                        size_start_stv = NaN;
                        starvation_entry_time = NaN;
                        first_divisiontime_stv = NaN;
                    end
                    break;
                end
                Gs = Gs +1;
            end
        end
        
        % Compute generations in recovery period (Gr)
        if bt +d < switchtime_2
            Gr = NaN;
        else
            Gr = 0;
            for  j = 1:length(thislineage)
                bt_start = timecourse(min(schnitzcells(thislineage(j)).frames));
                if (bt_start > switchtime_2)
                    Gr = Gr + 1;
                end
            end
        end
      
        % Estimating the size increase rate in each time zone by tracing
        % back to the lineage history.
        
        % Computing the size increase rate using the user-defined timezone
        % information for a given set of lineage info (thislineage).
        % Here, we only focus on the lineage history that already existed
        % before the recovery period
        if (bt <switchtime_2)
            [mdl_dt_s, ~] = compute_parameter_timezone(schnitzcells, thislineage, timecourse, timezone, ...
                delta_fr_stv, delta_t_stv, micronsperpixel);
        else
            mdl_dt_s = NaN(1,length(timezone)-1);
        end
        
    end

    %% Counting the number of offsprings before entering recovery period
    function countoffsprings_before_recovery(daughterone)
        last_time = data.birth_time(data.Schnitz == daughterone) + data.duration(data.Schnitz == daughterone);
        if last_time<=  switchtime_2
            if  ~isnan(data.daughter1(data.Schnitz == daughterone))
                countoffsprings_before_recovery(data.daughter1(data.Schnitz == daughterone));
                countoffsprings_before_recovery(data.daughter2(data.Schnitz == daughterone));
            end
        else
            offsprings = offsprings+1;
            
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
end