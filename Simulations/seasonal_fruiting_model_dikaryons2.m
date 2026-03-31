<--This code was used to generate figures {img:samplepath}, {img:timedist}, {img:timedist2} and {img:avgfreqstationary}-->
%% This models allele frequencies forward in time under the dikaryon fungal model with discrete fruiting seasons. 

function seasonal_fruiting_model_dikaryons2

%%PARAMTERS
alpha = 0.2; %rate of turnover in dikaryon population
s = 0.8; %rate of turnover in monokaryon population
F = 9;     %number of fruting bodies
M = 100;   %number of monokaryons
D = 100;   %number of dikaryons

time = 100;     %number of seasons
sim = 500;     %number of simulations
n_trait = 2;    %number of alleles
bone = ones(n_trait,1); %used in calculations

initY  = [0.3; 0.7];            %Initial data for the monokaryons 
initW = [0,0.5;0.5,0];  %Initial data for the dikaryons


%%INITIALISE VECTORS FOR STORING DATA
% tplot = [3, 10, 20, 50, 100, 400];  %Used when plotting the distribution at various times. Records the times at which to plot

Yavg = zeros(1,time);      % Mean allele frequency in monokaryons over all simulations
Wavg = zeros(1,time);       % Mean allele frequency in dikaryons over all simulations

Ytemp = zeros(1,1,time);      % Temporarily stores data to calculate Yavg and change dimensions
Wtemp = zeros(1,1,time);    % Temporarily stores data to calculate W11avg and change dimensions




%%BEGIN SIMULATION

for ii = 1:sim

     
[Y,W,Z] = Fmodel(s, alpha, F, M, D, time, initY, initW, bone, n_trait);



%For plotting sample paths:
% Y1(:,ii) = Y(1,1,:);
% dataY= Y1;

%For plotting mean allele values:
Ytemp = Ytemp + Y(1,1,:) ;        % Sum all monokaryon A allele frequencies pointwise in time
Wtemp = Wtemp + 2*W(1,1,:) + W(1,2,:) + W(2,1,:);       % Sum all dikaryon A allele frequencies pointwise in time
end

%For plotting mean allele values:
Yavg = Ytemp(1,1,:)./sim;      % Normalise Yavg to obtain an average
Wavg = 0.5*Wtemp./sim;  % Normalise W11avg to obtain an average




%% CALCULATE STATISTICS

%%FIGURES

% %%plot mon allele freqs through time
% figure
% plot(1:time, dataY)
% %title(['Frequency of allele A in monokaryons through time'])
% xlim([0 time])
% xlabel('Number of fruiting seasons')
% ylabel('Allele frequency')


% %% Plot distribution at various times
% figure
% for jj = 1:6 
%     subplot(2,3,jj)
%     histogram(Y1(tplot(jj),:),20, 'Normalization','probability')
%     xlim([0 1])
%     ylim([0 0.5])
% title(['time =' num2str(tplot(jj))])
% xlabel('Proportion of type A')
%  ylabel('Normalised frequency')
% end


%% Plot averages through time
figure
plot(1:time, Yavg(1,:))
hold on
plot(1:time, Wavg(1,:))
xlim([0 time])
ylim([0.3 0.6])
xlabel('Number of fruiting seasons')
ylabel('Average frequency')
legend('Monokaryons','Dikaryons')


end










function output = multidistmat(N,Z,ntrait)
%% This function takes a matrix-valued multinomial sample by converting the matrix of probabilities to a vector, and the vecto output back to a matrix
Z_temp = reshape(Z,[ntrait^2,1]);
out_temp = mnrnd(N, Z_temp);
output = reshape(out_temp, [ntrait,ntrait]);
   
end

function [Y,W,Z] = Fmodel(s, alpha, F, M, D, time, initY, initW, bone, n_trait) 
%% This function applies the fungal model of forward-time allele frequencies

Y = zeros(n_trait, 1, time);        %stores allele frequencies in monokaryons for each time in each simulation 
W = zeros(n_trait, n_trait, time);  %allele frequencies in dikaryons
Z = zeros(n_trait, n_trait, time);  %stores allele counts observed in sporocarps for each time in each simulation

Y(:,:,1) = initY;                           %apply initial data
W(:,:,1) = initW;                           %apply initial data
Z(:,:,1) = multidistmat(F,initW, n_trait);  %simulate initial sporocarps

for ii = 2:time
    
Q = (1-s)*Y(:,1,ii-1) + 0.5*s*(Z(:,:,ii-1)*bone+(bone'*Z(:,:,ii-1))')/F;    %Vector of monokaryon allele probabilities
Y(:,1,ii) = mnrnd(M, Q)/M;                                                  %Vector of realised monokaryon allele frequencies
P = (1-alpha)*W(:,:,ii-1) + alpha*Y(:,1,ii)*Y(:,1,ii)';                     %Matrix of dikaryon genotype probabilities
W(:,:,ii) = multidistmat(D, P, n_trait)/D;                                  %Matrix of realised dikaryon genotype frequencies
Z(:,:,ii) = multidistmat(F, W(:,:,ii), n_trait);                            %Matrix of sporocarp genotype counts

end
end
