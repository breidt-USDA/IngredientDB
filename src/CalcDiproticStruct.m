function res = CalcDiproticStruct(ConcAcid,pH,pKo1,pKo2,I)
%Fred Breidt
%Version 06102014
%CalcTriProtic
%Returns a "struct" with the pH, anion concentrations (Molar) and 
%the calculated total acid (as a check, it should equal ConcAcid) 
%See page 99 and page 168 of Butler and Cogley ("Ionic equilibrium", 
%1998, Wiley and Sons, Inc., NY). 
%
%Arguments:
%   Conc    Total concentration of the acid
%   pKo1    The pKo in pure water at 25C for the 1st pK (uncharged and -1)
%   pKo2    The pKo in pure water at 25C for the 2nd pK (-1 and -2)
%   I       Ionic strength of the solution 
%   adjC    Negative for cation, positive for anion
%           NOTE: if both cation and anion present, add together...
%The function adjusts pKa values for ionic strength using the 
%   pKa = adjpka(pKo, I) function
%   NOTE: Temp is igonored as effect is very small...
%Returns:
%   Struct with:
%   pH (pH value, Anion1(concentration of -1 charge anion form), 
%   Anion2(-2 charge), ProtA (protonated acid conc), 
%   and Total, the total acid concentration (sum of anions + prot)
%   NOTE, "res.Total" should equal the "Conc" parameter argument...
%   This serves as a great check that everything is working correctly. 
%NOTE: 
%   See CalcpH_AB.m for details with monoprotic acids...
%   See page 176 Butler and Cogley for diprotic acids...
%************************************************************************
H = 10^(-pH);
Kw = 1e-14;
pKa1 = adjpka(pKo1, I);          
pKa2 = adjpka(pKo2, I);
Ka1 = 10^(-pKa1);
Ka2 = 10^(-pKa2);
%An1 and An2 are the terms for each anion form...
%do this for each pH (i.e. H value)
denom = (H.^2 + (Ka1.*H) + (Ka1*Ka2));
An0 = ConcAcid * (H^2)/denom;
An1 = ConcAcid * (Ka1*H)/denom;
An2 = ConcAcid * (Ka1*Ka2)/denom;
%post results
res.Prot = An0;
res.Anion1 = An1; 
res.Anion2 = An2;
res.sum = An0+An1+An2;
end

function pKa = adjpka(pKo, IS)
%FUNCTION pKa = adjpka(pKo, IS)
%Returns the modified pKa of a monoprotic acid based on the
%ionic strength and a given temperature in degrees C. 
%See page 97 of Butler and Cogley ("Ionic equilibrium", 1998, 
%Wiley and Sons, Inc., NY). 
%
%Arguments:
%   pKo     The pKo in pure water at 25C or a vector of pKo's
%   IS      Ionic strength of the solution
%NOTE: we found that temperature does not significantly affect pKa
%by printing IS, temperature, and predicted pKa graph, so temperature 
%was not included.
%-------------------------------------------------------------------------          
b = 0.3;                            %from p97 of Butler and Cogley (B & C)
A = 0.51;                           %from p45 of B & C
SqrIS = sqrt(IS);                   %Square root of IS
temp = (SqrIS/(1 + SqrIS))-b*IS;    %equation 9 p97 B & C
pKa = pKo - (2*A*temp);             %equation 9 p97 B & C
end