function R51PolymerModelSigma_2

close all

number_of_solution_points = 50;

k_self_on = 10;  %nM^-1 * s^-1
k_self_off = 1; %s^-1
k_ssDNA_on = 10;   %nM^-1 * s^-1
k_ssDNA_off = 1;  %s^-1
k_dsDNA_on = 9.5;   %nM^-1 * s^-1
k_dsDNA_off = 1;  %s^-1

Rad51_start = 100; %nM (Concentration) (Estimate taken from Reuters et al. 2014)

number_of_hg_bps = 3.3e9;
number_of_euchromatic_bps = number_of_hg_bps * 0.5; %assume half of all DNA is accessible (euchromatin)
potential_Rad51_dsDNA_binding_sites = 0.8 * number_of_euchromatic_bps / 3; %assume 80% of DNA is in dsDNA form at any specific time
potential_Rad51_ssDNA_binding_sites = 0.2 * number_of_euchromatic_bps / 3; %assume 20% of DNA is in ssDNA form at any specific time
volume_of_cell = 2e-12; %L (Volume)
volume_of_nucleus = 0.1 * volume_of_cell; %L (Volume)
Avogadro_number = 6.023e23;

Rad51_ssDNA_binding_sites_concentration = ((potential_Rad51_ssDNA_binding_sites / Avogadro_number) / volume_of_nucleus) * 1e9; %nM - Estimate of concentration of ssDNA Rad51 binding sites (basepair triplets)...in average human cell
Rad51_dsDNA_binding_sites_concentration = ((potential_Rad51_dsDNA_binding_sites / Avogadro_number) / volume_of_nucleus) * 1e9; %nM - Estimate of concentration of dsDNA Rad51 binding sites (basepair triplets)...in average human cell

%%%%%%%%%%%%%%%%% Free Parameter %%%%%%%%%%%%%%%%
sigma = logspace(0, 10, number_of_solution_points); %local concentration effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_interval = 1e7;
number_of_species = 29;
number_of_observables = 20;

steadystate_values = zeros(number_of_observables, number_of_solution_points);
initial_conditions = zeros(1, number_of_species);
initial_conditions(1) = Rad51_start;

%matlabpool('open', 4);

for i = 1:number_of_solution_points
    parameters = [k_self_on, k_self_off, k_ssDNA_on, k_ssDNA_off, k_dsDNA_on, k_dsDNA_off, Rad51_start, Rad51_ssDNA_binding_sites_concentration, Rad51_dsDNA_binding_sites_concentration, sigma(i)];
    
    [t, y] = ode15s(@(time, species) R51Polymer2_Deriv(time, species, parameters), [0, time_interval], initial_conditions);
    steadystate_values(:, i) = R51Polymer2_Obs(y(end, :));
end

%matlabpool close
Colours = {'k','b','r','g','m',[.5 .6 .7],[.8 .2 .6]};

for i = 1:5
    subplot(2,3,1);
    plot(log10(sigma), Avogadro_number * (steadystate_values(i,:)/1e9), 'color', Colours{i}); %multiply steady state value by Avogadro's number and nucleus volume to get amount of R51 molecules in a specific multimeric state
    xlabel('log10(sigma)')
    out = ['Number of R51 Molecules'];
    ylabel(out)
    title('All Rad51 Multimers')
    hold on
end

legend('Monomers', 'Dimers', 'Trimers', 'Tetramers', 'Pentamers')

for i = 1:5
    subplot(2,3,2);
    plot(log10(sigma), i * steadystate_values(i + 10,:), 'color', Colours{i}); %multiply steady state value by i to get amount of R51 molecules in a specific multimeric state on ssDNA
    xlabel('log10(sigma)')
    out = ['Concentration (nM) * Oligomeric State'];
    ylabel(out)
    title('Rad51 Multimers on ssDNA Multiplied by Multimer Length')
    hold on
end

legend('Monomers', 'Dimers', 'Trimers', 'Tetramers', 'Pentamers')

for i = 6:10
    subplot(2,3,4);
    plot(log10(sigma), steadystate_values(i,:), 'color', Colours{i - 5}); %multiply steady state value by i to get amount of R51 molecules in a specific multimeric state
    xlabel('log10(sigma)')
    out = ['Concentration (nM)'];
    ylabel(out)
    title('Rad51 Multimers in Solution')
    hold on
end

legend('Monomers', 'Dimers', 'Trimers', 'Tetramers', 'Pentamers')

for i = 11:15
    subplot(2,3,5);
    plot(log10(sigma), steadystate_values(i,:), 'color', Colours{i - 10}); %multiply steady state value by i to get amount of R51 molecules in a specific multimeric state
    xlabel('log10(sigma)')
    out = ['Concentration (nM)'];
    ylabel(out)
    title('Rad51 Multimers on ssDNA')
    hold on
end

legend('Monomers', 'Dimers', 'Trimers', 'Tetramers', 'Pentamers')

for i = 16:20
    subplot(2,3,6);
    plot(log10(sigma), steadystate_values(i,:), 'color', Colours{i - 15}); %multiply steady state value by i to get amount of R51 molecules in a specific multimeric state
    xlabel('log10(sigma)')
    out = ['Concentration (nM)'];
    ylabel(out)
    title('Rad51 Multimers on dsDNA')
    hold on
end

legend('Monomers', 'Dimers', 'Trimers', 'Tetramers', 'Pentamers')


