% =========================================================================
% MATLAB SCRIPT: Run Optimization Problem in a Loop and Save Data
% =========================================================================
% This script performs a series of simulations by repeatedly solving a
% convex optimization problem. It generates random inputs for each
% iteration, calls a core CVX function to find the solution, and
% systematically stores all results. Finally, it saves the collected data
% in both .mat (for full fidelity) and .csv (for easy access) formats.
% =========================================================================

% --- Workspace Initialization ---
clear;         % Clear all variables from the workspace
clc;           % Clear the command window
close all;     % Close all open figure windows

% =========================================================================
% 1. Parameter Setup
% =========================================================================
N = 2;                  % Define the parameter N, related to system size.
num_iterations = 1000;  % Set the total number of iterations for the simulation.

% Calculate the length of the vector 'sig', which is 2^N.
sig_length = 2^N;

% =========================================================================
% 2. Pre-allocation for Results Storage
% =========================================================================
% Using a 'table' is a good practice for storing mixed data types clearly.
% Pre-allocating memory before a loop significantly improves performance in MATLAB.
results_table = table(...
    'Size', [num_iterations, 8], ... % Define table dimensions
    'VariableTypes', {               % Define data type for each column
        'double', ... % Iteration: The current loop number
        'double', ... % a: The random expectation value
        'cell',   ... % Input_sig: The input probability vector
        'double', ... % Cr_min: The minimized objective value from CVX
        'cell',   ... % O: The name of the selected Pauli operator
        'cell',   ... % Output_sig: The output probability vector from CVX
        'cell',   ... % pobj_data: History of the primal objective value
        'cell'      % dobj_data: History of the dual objective value
    }, ...
    'VariableNames', {               % Assign names to each column
        'Iteration', 'a', 'Input_sig', 'Cr_min', ...
        'O', 'Output_sig', 'pobj_data', 'dobj_data'
    } ...
);

% --- Script Start Notification ---
fprintf('Script started. Total iterations: %d.\n', num_iterations);
disp('--------------------------------------------------');
fprintf('Generating %d Pauli operators for N=%d...\n', 4^N, N);
pauli_operators = generate_paulis_with_names(N); % Generate Pauli operators and their names
num_paulis = length(pauli_operators);
fprintf('Generation complete.\n');

% =========================================================================
% 3. Main Simulation Loop
% =========================================================================
for i = 1:num_iterations
    
    % --- Step 3.1: Generate Random Inputs for the Current Iteration ---
    % Generate a random 4x4 density matrix.
    Rdrho = RandomDensityMatrix(4);

    % Randomly select one Pauli operator from the generated set.
    rand_idx = randi(num_paulis);
    selected_op_struct = pauli_operators(rand_idx);

    % Extract the matrix and its corresponding name from the struct.
    random_op_matrix = selected_op_struct.matrix;
    selected_op_name = selected_op_struct.name;
    
    % Calculate the expectation value of the operator on the density matrix.
    % This value is 'a'. real() is used to remove negligible imaginary parts
    % that can arise from floating-point inaccuracies.
    a = real(trace(random_op_matrix * Rdrho));
    b = a; % In this problem setup, b is equal to a.
    
    % Construct the operator matrix 'O' for the CVX function.
    O = ['II'; selected_op_name];
    
    % The input signal 'sig_input' is the diagonal of the density matrix.
    sig_input = diag(Rdrho);

    % --- Step 3.2: Prepare Inputs for the Optimization Function ---
    a_new = [1; a]; % Prepend 1 to 'a' as required by the function.
    b_new = [1; b]; % Prepend 1 to 'b' as required by the function.
    
    % Display progress in the command window.
    fprintf('Running iteration %d / %d... a = %.4f\n', i, num_iterations, a);
    
    % --- Step 3.3: Call the Core Optimization Function ---
    % This function solves the convex optimization problem.
    [Cr_min, sig_output, pobj_data, dobj_data] = ...
        CVX_Cr_with_sig_his(O, a_new, b_new, sig_input);
        
    % --- Step 3.4: Store All Results in the Pre-allocated Table ---
    % Use standard indexing for numeric types and curly braces {} for cell arrays.
    results_table.Iteration(i) = i;
    results_table.a(i) = a;
    results_table.Input_sig{i} = sig_input;
    results_table.Cr_min(i) = Cr_min;
    results_table.O{i} = selected_op_name;
    results_table.Output_sig{i} = sig_output;
    results_table.pobj_data{i} = pobj_data;
    results_table.dobj_data{i} = dobj_data;
    
end

% =========================================================================
% 4. Completion and Data Saving
% =========================================================================
disp('--------------------------------------------------');
fprintf('All %d iterations have been completed.\n', num_iterations);

% --- Step 4.1: Preview the Full Results Table ---
disp('Preview of complete results (MATLAB Table):');
disp(head(results_table, 10)); % Display the first 10 rows.

% --- Step 4.2: Prepare and Save a Numeric-Only CSV File ---
% This section converts the table to a "flat" format suitable for CSV.

% Extract only the scalar data columns into a new table.
scalar_results_table = results_table(:, {'Iteration', 'a', 'Cr_min', 'O'});

try
    % Convert the cell array of 'Input_sig' vectors into a single numeric matrix.
    % This syntax [cell{:}]' is an efficient way to concatenate vectors.
    % This will only work if all vectors in the cell array have the same length.
    sig_matrix = [results_table.Input_sig{:}]';
    
    % Get the number of elements in the signal vector to create column names.
    num_sig_elements = size(sig_matrix, 2);
    
    % Create new column names for the signal data (e.g., "sig_1", "sig_2", ...).
    sig_col_names = "sig_" + string(1:num_sig_elements);
    
    % Convert the numeric matrix into a table with the new column names.
    sig_table = array2table(sig_matrix, 'VariableNames', sig_col_names);
    
    % Concatenate the scalar data table and the new signal data table.
    final_csv_table = [scalar_results_table, sig_table];
    
    % Display a preview of the final, flattened table.
    disp('Preview of flattened, numeric-only table for CSV export:');
    disp(head(final_csv_table, 5));
    
    % Write the flattened table to a CSV file.
    csv_filename = 'simulation_results_numeric.csv';
    writetable(final_csv_table, csv_filename);
    fprintf('Flattened numeric results successfully saved to: %s\n', csv_filename);
    
catch ME
    % Error handling in case the cell-to-matrix conversion fails.
    fprintf('Error converting cell to double! Check if vector lengths are consistent.\n');
    fprintf('Error message: %s\n', ME.message);
    disp('Saving only scalar data as a fallback.');
    writetable(scalar_results_table, 'simulation_results_scalar_only.csv');
end

% --- Step 4.3: Save the Complete Dataset as a .mat File ---
% The .mat file preserves the original data structure, including cell arrays,
% which is crucial for a complete and accurate archive of the simulation.
mat_filename = 'simulation_results_complete.mat';
save(mat_filename, 'results_table');
fprintf('All complete results (including cell data) successfully saved to: %s\n', mat_filename);
disp('--------------------------------------------------');
