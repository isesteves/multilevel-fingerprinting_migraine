function [Idiff, Iself, Iothers] = diff_identifiability(corr_matrix, template)
% [Idiff, Iself, Iothers] = diff_identifiability(corr_matrix, template)
%
% Calculates differential identifiability based on correlation matrix and template.
%
% INPUTS:
%   - corr_matrix: Square correlation matrix representing pairwise correlations.
%   - template: Binary template indicating self and others (1 for self, 0 for others).
%
% OUTPUTS:
%   - Idiff: Differential identifiability, a measure of the difference between self and others.
%   - Iself: Mean correlation among elements marked as self in the template.
%   - Iothers: Mean correlation among elements marked as others in the template.
%
% USAGE:
%   [Idiff, Iself, Iothers] = diff_identifiability(corr_matrix, template)
%
% DESCRIPTION:
%   The function computes differential identifiability, which is a measure
%   of the difference in mean correlations between elements marked as self
%   and elements marked as others in the binary template. The template is
%   used to categorize elements in the correlation matrix into self and others.
%
% EXAMPLE:
%   % Generate a random correlation matrix and template
%   corr_matrix = randn(5, 5);
%   template = randi([0 1], 5, 1);
%
%   % Calculate differential identifiability
%   [Idiff, Iself, Iothers] = diff_identifiability(corr_matrix, template);
%
%   % Display the results
%   disp(['Differential Identifiability: ', num2str(Idiff)]);
%   disp(['Mean correlation among self: ', num2str(Iself)]);
%   disp(['Mean correlation among others: ', num2str(Iothers)]);
%

    % Calculate mean correlation among elements marked as self
    Iself = mean(corr_matrix(template == 1));
    
    % Calculate mean correlation among elements marked as others
    Iothers = mean(corr_matrix(template == 0));
    
    % Calculate differential identifiability
    Idiff = 100 * (Iself - Iothers);
end
