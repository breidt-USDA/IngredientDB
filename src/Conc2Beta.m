function BCmat = Conc2Beta( ABmat )
%Takes a concentration (in moles/liter) and a pH, and returns a BC value 
%  for that concentration at the given pH.
%DEPENDENCY: BetaModel_AB
%NOTE: this conversion does not consider salt concentration shifts in pK,
%  it just returns a beta value for a given pKa
%=========================================================================

%get number of rows in ABmat, set up return matrix
[rows, ~] = size(ABmat);            %rows in ABmat
NaClpercent = 0;                    %don't alter pK due to salt
temp = zeros(rows,2);               %return matrix

%for each row, calculate the beta value for the corresponding pH
for i=1:rows
   temp(i,:) = BetaModel_AB(ABmat(i,:),NaClpercent, ABmat(i,2));
end

%BetaModel returns pH, BC: we want BC, pH
BCmat(:,1) = temp(:,2);             %set beta value
BCmat(:,2) = temp(:,1);             %set pK value
end

