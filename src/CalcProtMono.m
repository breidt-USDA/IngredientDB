function ProtAcid = CalcProtMono( ISval, pHval, pKval, Conc_mM )
%Assumes 25C, monoprotic acid. calcs protonated acid for a given pH, pK, 
% ionic str. and acid concentration.

newpK = AdjustpKaMonoprotic(pKval, ISval, 25);

ProtAcid = Conc_mM/(1+10^(pHval - newpK));


end

