function qrs = phase2qrs(phase,debug)
% this function is used for finding the qrs locations using the phase generated 
% with the ecg model. A qrs is defined as a zero crossing.
%
% inputs
%   phase:  phase [sample]
%   debug:  enter debug? (default: 0) [bool]
%
% output
%   qrs:    qrs position [sample]
%
% fecgsyn toolbox, version 1.0, July 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 31-07-2014
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

% == check inputs
if nargin<2; debug=0; end;

% == zero phase crossing
sy = sign(phase);
flag_cross = (diff(sy)==2)|(sy(1:end-1)==0); % only from negative to positive (and not positive to negative) and zeros
qrs = find(flag_cross);

% == debug
if debug
   FONT_SIZE = 15;
   plot(phase,'LineWidth',3);
   hold on, plot(qrs,phase(qrs),'+r','LineWidth',2);
   legend('phase','QRS');
   xlabel('Time [sec]'); ylabel('FHR [bpm]')
   set(gca,'FontSize',FONT_SIZE);
   set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
end

end