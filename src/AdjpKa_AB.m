function ABmat = AdjpKa_AB(inABmat, IS)
%From a matrix of conc, pKa adjust pKa's for ionic strength (IS)
% NOTE: Typically this is NaCl, but IS may be calculated for any salts
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

temp = 25;                          %temp in C, by default 25C
[rows,~] = size(inABmat);           %get rows
ABmat = inABmat;                    %default output with initial pKs

%calculate and assign the new pKa values for each row 
for i=1:rows
    ABmat(i,2) = AdjustpKaMonoprotic(ABmat(i,2),IS,temp);
end


end

