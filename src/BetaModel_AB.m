function BCcurve = BetaModel_AB(ABmatrix,NaClpercent,pHvec)
%Parameters
%   ABmatrix = Nx2 matrix with columns of Conc (Molar), pKa
%   NaClpercent = Percent NaCl in solution (needed for ionic str. calc)
%   pHvec = Nx1 col vector of pH values at which to determine B val 
%NOTE:  assumes equal number of conc. values and pKas
%RETURNS
%   BuffCapCurve = Nx2 array (pH, BCvalue)
%Model: Butler and Cogely, 1998. Ionic Equilibrium, Solubility and pH 
%   Calculations. John Wiley and Sons, Inc. NY. Equation 82 on page 134.
%SAMPLE CALL: acetic and lactic at 0.02 M with pKas = 3.86
%   >> BC = BetaModel([0.05 4.76], 2, [2:0.1:12]');
%   >> plot(BC(:,1),BC(:,2),'or')
%NOTE: NaCl can be 0 if not in a salt solution for a close
%   approximation of the answer.
%Const can have values:
%const = 2.302585;      ln(10) constant
%--------------------------------------------------------------------------

%Conversions and assignments
Conc = ABmatrix(:,1);                       %concentraitons for buffers
pKa = ABmatrix(:,2);                        %pKa values for buffers
const = 2.302585;                           %Natural log of 10
Temp = 25;                                  %assume 25C for temperature
IonicStr = NaClpercent/5.84;                %convert NaCl to IS
Kw = 10^-(AdjustpKaMonoprotic(14,IonicStr,Temp)); %water equilib. constant
pKa =AdjustpKaMonoprotic(pKa,IonicStr,Temp); %adjust pKa's
Ka =  10.^(-1*pKa);                         %convert pKa to Ka
Nbuffers = max(size(Conc));                 %number of buffers
H = 10.^(-1*pHvec);                         %convert pH vector to [H+]
NpHvals = max(size(pHvec));                 %number of pH values in vector
OH = Kw./H;                                 %calc OH from H and Kw
bufvecs = zeros(NpHvals,Nbuffers);          %matrix to hold buffer deriv
BCcurve = zeros(NpHvals,2);                 %buffer capacity curve (Nx2)

%sum deriv of acid base terms (from Butler and Cogley, eq. 82 on p.134)
for j = 1:Nbuffers                          %for each buffer
    bufvecs(:,j) = ((Conc(j)*Ka(j))*H)./((H + Ka(j)).^2); %calc derivative
end
bufsum = sum(bufvecs,2);                    %sum rows

%return matrix (Nx2) pH values, beta value
BCcurve(:,1) = pHvec;                  %assign pH values
BCcurve(:,2) = const*(bufsum + OH + H);%eq. 82, p.134 of B&C

%SUBFUNCTION==============================================================
function pKa = AdjustpKaMonoprotic(pKo, I, TempC)
%FUNCTION pKa = AdjustpKaMonoprotic(pKo, I, TempC)
%Returns the modified pKa of a monoprotic acid based on the
%ionic strength and a given temperature in degrees C. 
%See page 97 of Butler and Cogley ("Ionic equilibrium", 1998, 
%Wiley and Sons, Inc., NY). 
%
%Arguments:
%   pKo     The pKo in pure water at 25C or a vector of pKo's
%   I       Ionic strength of the solution 
%   TempC   The temperature in degrees C
%   Z       Vector of ionic charges
%   b       Value of adjustment parameter for Davies equation = 0.3
%-------------------------------------------------------------------------          
b = 0.3;                %adjustment factor from p97 of Butler and Cogley
Epsilon = 78.3808;                          %dialectric constant of water
degK = TempC + 273;                         %calc temp in K
A = 1.825 * 10^6 * (Epsilon * degK)^(-3/2); %from p45 of Butler and Cogley
temp = (sqrt(I)/(1 + sqrt(I)))-(b*I);       %equation 9 p97 B & C
pKa = pKo - (2*A*temp);                     %equation 9 p97 B & C
end

end