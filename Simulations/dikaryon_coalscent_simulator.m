function coal_times = dikaryon_coalescent_simulator(coalescent_info_init, params, nsim)

%coalescent_info_init = [1,1,0];    %The initial states of each nucleus. This has coding
    % 0 = monokaryon
    % 1 = dikaryon
    % 2 = sporocarp
    % 501, 502, 503, ... = paired dikaryons, where the digits after 5 record the index position of the paired nucleus
    % 901, 902, ... = coalesced, where the digits after 9 record the position of the representative lineage which remains
 % The final entry is the time at which the configuration is observed (initially 0)
 % After every coalescence another row is added to this array, recording the new configuration
  
n = length(coalescent_info_init)-1;  % Initial number of lineages
%params = [0.4,0.5,0.7,20]; % betam, betad, s, F. Assume M, D infinite


coal_times = zeros(nsim,1);

for ii = 1:nsim
    coalescent_info = dikaryon_coalescent_sim(coalescent_info_init, n, params);
    coal_times(ii) = coalescent_info(2,3);
end


end




function coalescent_info = dikaryon_coalescent_sim(init, n, params)
% simulates the coalescence of n nuclei with an initial configuration coalescent_info_init under a model with parameters encapsulated by params
% params = (beta_m, beta_d, s, F)

coalescent_info = init; %Initialise the output with the initial data
states = init(1:n);  % Initial states of the nuclei
time = 0;




lineages = n;

    while lineages > 1
        %keep doing stuff until all lineages have coalesced
        time = time +1; %update the time/seasons

        % Apply the three functions which descibe how muclei move through states
        states = mating_func(states, params);
        states = spore_func(states, params(3));
        [states, flag] = sporocarp_func(states, params(4));
        
        if flag >0
            lineages = lineages - flag; % Flag records the number of lineages which coalesce
            newdata = [states, time];   % Add the time of coalescence to the state info
            coalescent_info = [coalescent_info; newdata];   % Add this new coalescence data to the output
        end


    end
end

% 0 = monokaryon
    % 1 = dikaryon
    % 2 = sporocarp
    % 10, 11, 12, ... = paired dikaryons, where both nuclei of the pair have the same value
    % 901, 902, ... = coalesced, where the digits after 9 record the position of the representative lineage which remains

    
function [states_out, flag] = sporocarp_func(states, F) %Yo this isn't doing the right thing
    %% Simulates where each nucleus in a sporocarp originated from
    % I apply this to each nucleus in turn
    index_nuclei = find(states(:) == 2);         % Indices of nuclei that are in sporocarps. 
    n_nuclei = length(index_nuclei);             % Number of nuclei that are in sporocarps. 
    flag = 0;
    sps = zeros(n_nuclei, 3);
    sps(:,3) = index_nuclei;
    states_out = states;
    
    for jj = 1:n_nuclei  % For each sporocarp nucleus
        sps(jj,1) = randi(F); % Assign it a numbered sporocarp of origin
        sps(jj,2) = randi(2); % Assign it a numbered nucleus of origin
    end

    for kk = 1:F  % For each sporocarp
        F_ind = find(sps(:,1) ==kk); %Make an array with the jj indices of all nuclei from sporocarp kk
        if length(F_ind) == 1 % No other nuclei came from this sporocarp
            states_out(sps(F_ind,3)) = 1; % comes from a dikaryon
        elseif length(F_ind) > 1 % At least one other nucleus is from the same sporocarp
            % Is there a paired dikaryon going on?
            if sum(sps(F_ind,2)==1)==0 || sum(sps(F_ind,2)==2)==0 %If only one nucleus is a parent
                %1st nucleus came from a dikaryon, all the others coalesce
                states_out(sps(F_ind(1),3)) = 1;     % The diakaryon
                for ll = 2:length(F_ind)
                    states_out(sps(F_ind(ll),3)) = 900+sps(F_ind(1),3);  % Coalesced, with a representative lineage
                    flag = flag+1;
                end
            else % Paired dikaryon
                %This is harder to track
                    %I want to get the indices of all the 1 nuclei and the
                    %indices of all the 2 nuclei, and repeat
                n1 = find(sps(F_ind,2)==1);
                n2 = find(sps(F_ind,2)==2);
                
                states_out(sps(F_ind(n1(1)),3)) = 500+ sps(F_ind(n2(1)),3);
                states_out(sps(F_ind(n2(1)),3)) = 500+ sps(F_ind(n1(1)),3);
                for mm = 2:length(n1)
                    states_out(sps(F_ind(n1(mm)),3)) = 900+sps(F_ind(n1(1)),3);  % Coalesced, with a representative lineage
                    flag = flag+1;
                end
                for nn = 2:length(n2)
                    states_out(sps(F_ind(n2(nn)),3)) = 900+sps(F_ind(n2(1)),3);  % Coalesced, with a representative lineage
                    flag = flag+1;
                end
            end
            
                
        end
        
    end

        
end



    
function states = spore_func(states, s)
    %% Simulates where each nucleus originated from over one influx of spores
    % I apply this to each nucleus in turn
    
    for ii = 1:length(states)        
        if states(ii) == 0 %it's a monkaryon
            if rand() < s
                states(ii) = 2; %came from a spore from a sporocarp
            end
        end
    end
end    
    
    
    
    

function states_out = mating_func(states, params)
    %% Simulates where each nucleus had originated from over a mating period
    % I apply this to each nucleus in turn
    betam = params(1);
    betad = params(2);
    alpha = betam+0.5*betad;
    n = length(states);
    states_out = states;
    
    skip = zeros(n,1);
    for ii = 1:n
        if skip == 0 %Skip any nuclei indexed 1 - this avoids dealing with paired nuclei twise
            
            if states(ii) == 1 %it's a dikaryon
                if rand() < alpha
                    states_out(ii) = 0; %came from a monokaryon by mating
                end
            elseif (states(ii) > 500) && (states(ii) < 600) %it's part of a dikaryon pair
                temprand = rand();
                if temprand < betam %it comes from mon-mon mating
                    states_out(states(ii)-500) = 0;
                    skip(states(ii)-500) = 1; %skip the other nucleus
                    states_out(ii) = 0;
                elseif temprand < betam+betad %it comes from di-mon mating
                    states_out(states(ii)-500) = 0;
                    skip(states(ii)-500) = 1; %skip the other nucleus
                    states_out(ii) = 1;
                end
            end
        end
    end
end

    
