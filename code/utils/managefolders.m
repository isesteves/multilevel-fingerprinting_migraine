function managefolders(foldername, method)
% managefolders(foldername, method)
%
% INPUT
% > foldername - (string) foldername, including path
% > method - (string) Options:
% - 'create' - checks if the folder exists, if not, it creates it

if strcmp(method, 'create')
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
end
