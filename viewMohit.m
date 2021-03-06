%% Quickly visualizes raw data to find channels of interest.
% will load all .mat files in a folder and generate plots for each. Figure
% numbers do not necessarily follow order in directory! Run with figDelete
% to get rid of bad files

Fs = 600e3;
files = dir(pwd);
f = ~[files.isdir];
subFolders = files(f);
fileNames = {subFolders.name}';

for i = 1:size(fileNames)
    currentFile = char(fileNames(i));
    V = load(currentFile);
    V1 = struct2cell(V);
    V = V1{1};
    V = double(V);
    x = (1:1:length(V))/Fs;
    figure
    plot(x,V);
    title(currentFile)
end