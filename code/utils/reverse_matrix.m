function reversed_matrix = reverse_matrix(original_matrix, max_val)

% Reverse the order of unique values
reversed_values = max_val:-1:1;

% Create a mapping using arrays
mapping = [0 reversed_values];

% Apply the mapping to reverse the order of integer values
reversed_matrix = arrayfun(@(x) mapping(x+1), original_matrix);
