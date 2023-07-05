function res = SCBCfit_area(BCcurve, order)
%Fit obs data with the sin-cos function, return the predicted fit
%  and an AB matrix with BC values and pKs for the peaks of BC found 
%  in the predicted data. 
%concentrations for BC model from sine-cosine fitted peaks
%PARAMETERS:
%  BCcurve = Nx2, observed BC data, col1 = pH, col2 = BC value
%  order = number of paired sine-cosine terms (harmonics)
%RETURNS
%  res = data structure
%  res.Obs = Nx2, original BC curve: pH, BC value
%  res.Pred = Nx2, predicted (SCBC) model data: pH, BC value 
%  res.SSE = sum sq. error term
%  res.recip_cond = matrix condition 
%  res.determinant = matrix determinant
%NOTE: mul = multiplier for X vals in sine-cosine model, set to pi/2 minus
%  a small amount. Trial and error found that 0.5 works pretty well. This
%     allows a bit less than one cycle of sine-cosine for pH (x) = 2..12
%*************************************************************************

%assignments
X = BCcurve(:,1);                %init X vals (Nx1) N= no. data points)
Y = BCcurve(:,2);                %init Y vals (Nx1)
mul = 0.5;                       %X val multiplier (for 1 cycle 2..12)                      
modelX = linspace(min(X),max(X),100)'; %pred X: 100 pts between 2..12\
% modelX = linspace(2,12,100)';    %actually use 2 to 12 for the model
terms = (2*order) + 1;           %b0 + 2 terms (for "degree" or "harmonic")
N = length(X);                   %number of X data points

%set up partial deriv matrix with summed sin, cos terms 
factormat = GetFactorMat(terms,order,N); %get matrix of X for each deriv
pmat = zeros(terms,terms);       %square matrix for partial derivatives
for l=1:terms                    %for each row of square matrix (L not 1)
   for m=1:terms                 %for each col of square matrix
      pmat(l,m) = sum(factormat(l,:).*factormat(m,:)); %sum sin and/or cos
   end  
end

%get Y vector
ymat = zeros(terms,1);           %col vector for Y vals * each term
for n=1:terms                    %for each term mult Y val * sin or cos
   ymat(n) = sum((factormat(n,:))'.*Y); %sum: SorC (col) * Y vals (col)
end

%matrix division to calc params for sine-cosine model
params = pmat\ymat;              %divide Y matrix by square param matrix

%calc model with final params and predicted X values
modelY = CalcSCModel(params,modelX); %calc predicted vals for each model X
SSE = CalcSC_SSE(params);        %use predicted parameters to calc SSE
area = CalcArea(params,min(X),max(X),mul);

%output
res.Obs = BCcurve;               %BC curve
res.Pred = [modelX modelY];      %predicted (SCBC) model data
res.params = params;             %parameter vector
res.area = area;                 %area under curve between xmin and xmax
res.SSE = SSE;                   %sum sq. error term
res.recip_cond = rcond(pmat);    %matrix condition
res.determinant = det(pmat);     %determinant

%SUBFUNCTIONS-------------------------------------------------------------

%FUNCTION: get a vector (termsxN) with 2*order + 1 number of sine - cosine 
% terms (with initial 1) and appropriate multiplier
function mat = GetFactorMat(terms,order,N)
   k = 2;                        %matrix row index, skip first row (all 1s)
   mat = zeros(terms,N);         %one row of X vals for each term 
   mat(1,:) = 1;                 %first row is all 1s 
   for j=1:order                 %for each set of terms (order)
      mat(k,:) = sin(j*mul*X);   %sine term (times order index and mult.)
      k = k+1;                   %advance row index
      mat(k,:) = cos(j*mul*X);   %cosine term (times order index and mult.)
      k = k + 1;                 %advance row index
   end
end

%FUNCTION: calc Y vals with a given param set and set of X vals w/sin-cos
function predY = CalcSCModel(currPrms,Xvals) %calc model for given coef, X
   loopvar = (length(currPrms)-1)/2;  %loop for each set of terms (but not Bo)
   b0 = currPrms(1);                %set Bo
   currPrms(1) = [];                %remove Bo param, leaves even no. of sets
   p = 1;                           %counter for param index
   sumterms = zeros(length(Xvals),1);  %variable for sum of sine-cosine terms
   for i = 1:loopvar                %for each set of sine-cosine terms
      sumterms = sumterms + currPrms(p)*(sin(i*mul*Xvals))...  
      + currPrms(p+1)*(cos(i*mul*Xvals));  %calc the vector for each term
      p = p + 2;                    %advance index to next set
   end 
   predY = b0 + sumterms;           %get pred Y values for each X
end

%FUNCTION: calc SSE with a given param set for sine-cosine model
function SSE = CalcSC_SSE(currPrms) %calc sum squared error for param set 
   predYvals = CalcSCModel(currPrms,X);  %uses obs X values for this
   SSE = sum((predYvals-Y).^2);  %sum squared errors, predicted minus Y      
end

function area = CalcArea(params,minX,maxX,alpha)
   loopvar = (length(params)-1)/2;  %loop for each set of terms (but not Bo)
   b0 = params(1);                  %set Bo
   params(1) = [];                  %remove Bo param, leaves even no. of sets
   p = 1;                           %counter for param index
   evalMax = 0;                     %to sum max X answer
   evalMin = 0;                     %to sum min X answer
   for i=1:loopvar                  %loop for each set of terms
      evalMax = evalMax + (params(p+1)/(i*alpha))*sin(i*alpha*maxX) - ...
         (params(p)/(i*alpha))*cos(i*alpha*maxX);
      evalMin = evalMin + (params(p+1)/(i*alpha))*sin(i*alpha*minX) - ...
         (params(p)/(i*alpha))*cos(i*alpha*minX);
      p = p+2;
   end
   evalMax = (b0*maxX) + evalMax;
   evalMin = (b0*minX) + evalMin;
   area = evalMax - evalMin;
end

end

