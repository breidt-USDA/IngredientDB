function pKnew = AdjpKa(pKold, IS)
%Take a matrix of conc, pKa and adjust pKa's for salt (IS)
%   
%FB Version 08022018
%PURPOSE 
%   Adjust pKas in matrix for Ionic strength (assumes 25C)
%NOTE:
%   Assumes all acids/bases in matrix are monoprotic.
%PARAMETERS
%   inABmat = matrix (Nx2) of:  conc (in moles/L), pKa values (pH units)
%       NOTE: only one acid or base (conc, pKa pair) is needed.
%   IS = scalar, equalt to ionic strength (IS units)
%DEPENDENCIES
%   AdjustpKaMonoprotic
%-------------------------------------------------------------------------

%params
temp = 25;                          %temp in C


%calc new pKas
pKnew = AdjustpKaMonoprotic(pKold,IS,temp);


end

