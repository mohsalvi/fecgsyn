%% Main script for FECG morphological analysis
% 
% This is the mains script for testing the morphological consistency of
% extracting the foetal signal using various methods. 
% Used extraction methods:
%  - ICA
%  - PCA
%  -piCA
% 
% Used morphological measures:
% - T/QRS ratio
% - ST segment
% - QT interval
% 
% 
% 
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 03-06-2014
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Input parameters
% saving path
if isunix
    path = ['~/foobar/fecgdata/' datestr(date,'yyyy.mm.dd') '/'];
else
    path = ['C:\foobar\fecgdata\' datestr(date,'yyyy.mm.dd') '\'];
end

generate = 1;   % boolean, data should be generated? 
                % If not, path should direct to data location

%% Data Generation
if generate
    mkdir(path)
    generate_data(path)  % generates set of unique simulated data for testing
else
    cd(path)
end

%% Extraction Methods
fls = dir('*.mat');     % looking for .mat (creating index)
fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);
for i = 1:length(fls)
    load(fls)
    script_ica_pca(data)
    
end
%% Morphological Analysis

