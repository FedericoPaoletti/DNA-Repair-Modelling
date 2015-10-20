close all

n_points_start = 100;
n_points_sigma = 100;
k_RAD51_sson = 0.01; %Assume K_D (protein-ssDNA interaction) = 1000 microM = 0.001 M Set on rate of RAD51 for ssDNA to low value! Assume low affinity for DNA!
k_RAD51_ssoff = 10;
k_BRCA2_sson = 0.01; %Assume K_D (protein-ssDNA interaction) = 1000 microM = 0.001 M
k_BRCA2_ssoff = 10;
k_RAD51_dson = 0.01; %Assume K_D (protein-dsDNA interaction) = 1000 microM = 0.001 M
k_RAD51_dsoff = 10;
k_BRCA2_dson = 0.01; %Assume K_D (protein-dsDNA interaction) = 1000 microM = 0.001 M
k_BRCA2_dsoff = 10;
kpp_on = 0.01;
kpp_off = 0.0001; % Assume K_D (protein-protein interaction) = 0.01 M
BRCA2_start = logspace(0, 6, n_points_start); %microM
RAD51_start = logspace(0, 6, n_points_start); %microM
ssDNA_start = 100; %microM 
dsDNA_start = 100; %microM 

sigma = logspace(0, 8, n_points_sigma);

n_time = 1e7;
steadystate_value = zeros(1, n_points_start);
options = odeset('AbsTol', 1e-9, 'RelTol', 1e-6);
matlabpool('open', 4)

parfor j = 1:n_points_sigma
    j
    for i = 1:n_points_start
        parameters = [k_RAD51_sson k_RAD51_ssoff k_BRCA2_sson k_BRCA2_ssoff k_RAD51_dson k_RAD51_dsoff k_BRCA2_dson k_BRCA2_dsoff kpp_on kpp_off ssDNA_start dsDNA_start BRCA2_start(i) RAD51_start(i) sigma(j)];
        [t, y] = ode15s(@(time,species) integ_dsDNA_simple_improved_mex(time, species, parameters), [0 n_time], [ssDNA_start dsDNA_start BRCA2_start(i) RAD51_start(i) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0], options);
        %'steady state y values', y(end,:)
        steadystate_value(:, i, j) = obs_dsDNA_simple_improved(y(end,:));
    end
end

matlabpool close
%size(steadystate_value(1, :, :))
Bound_RAD51 = squeeze(steadystate_value(1, :, :));
normalised_RAD51 = zeros(length(RAD51_start));
Bound_RAD51
for i = 1:length(RAD51_start)
    normalised_RAD51(:, i) = Bound_RAD51(:, i)/(max(Bound_RAD51(:, i)));
end

surf(log10(sigma), log10(RAD51_start), normalised_RAD51); 
xlabel('log_{10}(sigma)')
ylabel('log_{10}(RAD51 start)')
zlabel('Normalised Bound RAD51')
figure
surf(log10(sigma), log10(RAD51_start), Bound_RAD51); %Plot all ternary complexes (1)
xlabel('log10(sigma)')
ylabel('log10(RAD51 start)')
zlabel('Bound RAD51')

%%%%%%%%HILL NUMBER DEFINITION ACCORDING TO GOLDBETER ET AL (1981) PNAS
%%%%%%%%PAPER
hill_number_2protein = zeros(100,1);
for i = 1:100
    
    [S_90_value, S_90_index] = min(abs(Bound_RAD51(:,i) - (max(Bound_RAD51(:,i))) * 0.9));
    [S_10_value, S_10_index] = min(abs(Bound_RAD51(:,i) - (max(Bound_RAD51(:,i))) * 0.1));   


    hill_number_2protein(i) = log(81)/(log(RAD51_start(S_90_index)/RAD51_start(S_10_index)));
end

hill_number_2protein

figure
plot(log10(sigma), hill_number_2protein, 'k')
xlabel('log_{10}(\sigma)')
ylabel('Hill Coefficient')
