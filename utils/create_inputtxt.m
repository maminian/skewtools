function create_inputtxt(filename,aratios,maxM,maxN,flowtype)
% Script for creating 'input.txt' with desired range of values.
%
% function create_inputtxt(filename,aratios,maxM,maxN,flowtype)
%
%
%


nAratios = length(aratios);

fid = fopen(filename,'w');

% Number of aspect ratios to use
fprintf(fid,'%i%s\n',nAratios,' ! Number of aspect ratios');

% List of aspect ratios
fprintf(fid,'%e, ',aratios);
fprintf(fid,'%s\n',' ! List of aspect ratios');

% Max sum index in M
fprintf(fid,'%i%s\n',maxM,' ! Max index in M');

% Max sum index in N
fprintf(fid,'%i%s\n',maxN,' ! Max index in N. Multiplied with M for index set sum.');

fprintf(fid,'%s',flowtype)

fclose(fid);

end
