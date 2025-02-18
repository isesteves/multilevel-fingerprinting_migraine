function data_out= conversion_vecmat(data_in, nr_areas, conversion)
% conversion_vecmat(data_in, nr_areas, conversion) 
%
% Converts vector into symmetric matrix or vice-versa
%
% INPUTS:
%   - data_in: Input data, either a vector or a 3D matrix.
%   - nr_areas: Number of areas (dimension of the square matrix).
%   - conversion: Conversion type, either 'vec2mat' or 'mat2vec'.
%
% OUTPUT:
%   - data_out: Output data, either a symmetric matrix or a vector.
%
% USAGE:
%   conversion_vecmat(data_in, nr_areas, 'vec2mat') converts a vector
%   'data_in' into a symmetric matrix 'data_out' with dimensions 'nr_areas x nr_areas'.
%
%   conversion_vecmat(data_in, nr_areas, 'mat2vec') converts a symmetric matrix
%   'data_in' into a vector 'data_out', where only the lower triangular part
%   excluding the main diagonal is considered.


if strcmp(conversion, 'vec2mat')
    
    nr_col = size(data_in, 2);
    data_out = zeros(nr_areas, nr_areas, nr_col);
    for n = 1:nr_col
        aux = tril(ones(nr_areas),-1);
        aux(aux > 0) = data_in(:,n);
        data_out(:,:,n) = (aux + aux')./(eye(nr_areas)+1);
    end  
    
elseif strcmp(conversion, 'mat2vec')  
%     tril_m  = tril(true(size(data_in)), -1);
%     data_out  = data_in(tril_m).';
    
     % Get the size of the input matrix
    [n, ~, m] = size(data_in);
    
    % Initialize the output vector
    data_out = zeros((n * (n - 1) / 2), m);
    
    % Loop over each slice along the third dimension
    for k = 1:m
        % Extract the k-th matrix
        currentMatrix = data_in(:,:,k);
        
        % Extract the lower triangular part excluding the main diagonal
        lowerTriangular = currentMatrix(tril(true(n, n), -1)).';
        
        % Assign the values to the corresponding column in the output vector
        data_out(:, k) = lowerTriangular;
    end
end


