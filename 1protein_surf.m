close all

n_points_kpp = 100;
n_points_kss = 100;
k_X_sson = 0.01; %Set on rate of RAD51 for ssDNA to low value! Assume low affinity for DNA!
k_X_ssoff = logspace(-7, 6, n_points_kss);
k_X_dson = 0.01; %Set on rate of RAD51 for ssDNA to low value! Assume low affinity for DNA!
k_X_dsoff = logspace(-7, 6, n_points_kpp);
X_start = 1; %microM
ssDNA_start = 10; %microM
dsDNA_start = 1000; %microM

n_time = 1e8;
steadystate_value = zeros(2, n_points_kpp);
options = odeset('AbsTol', 1e-9, 'RelTol', 1e-6);

matlabpool ('open', 4);

tic
parfor j = 1:n_points_kss
    j
    for i = 1:n_points_kpp
        parameters = [k_X_sson k_X_ssoff(j) k_X_dson k_X_dsoff(i) X_start ssDNA_start dsDNA_start];
        initial_conditions = zeros(1, 5);
        initial_conditions(1) = ssDNA_start;
        initial_conditions(2) = dsDNA_start;
        initial_conditions(3) = X_start;
        [t, y] = ode15s(@(time,species) integ_1protein_mex(time, species, parameters), [0 n_time], initial_conditions, options);
        %'steady state y values', y(end,:)
        steadystate_value(:, i, j) = obs_1protein(y(end,:));
    end
end
toc

matlabpool close;

%size(steadystate_value(1, :, :))
ssDNA_complexes = squeeze(steadystate_value(1, :, :));
dsDNA_complexes = squeeze(steadystate_value(2, :, :));

surf(log10(k_X_ssoff/k_X_sson), log10(k_X_dsoff/k_X_dson), ssDNA_complexes); %Plot all ternary complexes (1)
xlabel('log10(ssDNAKd)')
ylabel('log10(dsDNAKd)')
zlabel('ssDNA Complex Concentration')

figure
surf(log10(k_X_ssoff/k_X_sson), log10(k_X_dsoff/k_X_dson), dsDNA_complexes);
xlabel('log10(ssDNAKd)')
ylabel('log10(dsDNAKd)')
zlabel('dsDNA Complex Concentration')
