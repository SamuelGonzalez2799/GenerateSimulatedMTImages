function varargout = getFluorCoords(numUnits,tagScheme,schemeValues)
%fluorNumPerUnit = getFluorsCoords(numTubulin,tagScheme,schemeValues)
%   Returns the number of fluorophores per unit as a vector with one
%   element per unit
%[fluorNumPerUnit,fluorIndexIntoUnits] = getFluorsCoords(numTubulin,tagScheme,schemeValues)
%   The additional output returns a vector with length equal to the total
%   number of fluorophores, each element containing the index of the unit
%   to which the fluorophore is bound. 
%INPUTS:
%   numUnits: the number of units to simulate fluorophore incorporation
%   tagScheme: use one of the following options (input as a string, 
%       not case sensitive):
%       "none" : no fluorophores are incorporated on any units
%       "single" : units either incorporate one fluorophore or no
%           fluorophores, with a fixed probability of incorporation
%       "uniform" : units incorporate 0-nMax fluorophores with a uniform
%       probability distribution
%       "poisson" : units incorporate fluorophores based on a poisson
%           distribution
%       "custom" : a custom distribution, defined by the schemeValues.
%       "mixed" : allows for mixing any number of different taging schemes,
%           each defined scheme to be mixed allows a mixing ratio to be 
%           defined
%   schemeValue: a value or cell array of parameters, differs per scheme:
%       "none" : no scheme values required
%       "single" : (1) probability of incorporation (a value between 0-1)
%       "uniform" : (1) maximum possible per unit (integer value > 0)
%       "poisson" : (1) lambda (positive real number, describes expected value  
%           and variance of a poisson distribution), (2)(Optional) A value
%           determining the limit of number per unit: either a value
%           greater than 1 to describe the largest number per unit, or a
%           value less than 1 to describe probability density value beyond
%           which to ignore. 
%       "custom" : (1) a vector of N values, where the vector corresponds to 
%           the relative probability of a unit having 0 - (N-1)
%           fluorophores associated.
%       "mixed" : an N x(M+2) cell array. N is the number of distributions to be
%           mixed. The first column should contain strings defining the
%           tagScheme. The subsequent column contains the relative mixing
%           proportion of that scheme should contain whatever values are
%           needed for the given tagging scheme. M will be the largest
%           number of values required for any of the tag schemes used. Any 
%           schemes requiring fewer than M values can be left empty.
%           An example would be:
%           {'none',0.75,[];
%            'single',0.15,0.5;
%            'poisson',0.10,4}
%           This would indicate that the final output should contain 75%
%           non-tagged, 15% having a 50% chance of incorporation, and 10%
%           incorporating based on a poisson distribution with lamda = 4
%     

nargoutchk(1,2);

schemeAndVarsNeeded = {'none',0;'single',1;'uniform',1;'poisson',1;'custom',1};
%Mixed is not included, as all are converted to a "mixed" format.

if(strcmpi(tagScheme,'mixed'))
   if( iscell(schemeValues))
       if (size(schemeValues,2) < 2)
           error('schemeValues must contain at least 2 columns, one for the tagScheme string, and one for the relative proportion')
       end
       numDists=size(schemeValues,1);
   else
       error('Scheme Values must be a cell array for the tag scheme "mixed"');
   end
else
    schemeValues = {tagScheme,1,schemeValues};
    numDists=1;
end

for n = 1:size(schemeValues,1)
    thisScheme = strcmpi(schemeValues(n,1),schemeAndVarsNeeded(:,1));
    numVarsNeeded=schemeAndVarsNeeded{thisScheme,2};
    if(~size(schemeValues,2)>= numVarsNeeded+2)
        error('The schemeValues does not contain enough columns for a scheme of type %s',schemeAndVarsNeeded{thisScheme,1});
    end
end

%figure out the relative proprtions and convert to a normalized value
totalProp = sum(cell2mat(schemeValues(:,2)));
schemeValues(:,2)=num2cell(cell2mat(schemeValues(:,2))/totalProp);

probDist(1) = 0;

for d = 1:numDists
    switch lower(schemeValues{d,1})
        case 'none'
            proportion = schemeValues{d,2};
            probDist(1)=probDist(1)+proportion;
            
        case 'single'
            proportion = schemeValues{d,2};
            prob = schemeValues{d,3};
            if(numel(probDist)<2)
               probDist(end+1:2) = 0; 
            end
            probDist(1:2) = probDist(1:2) + ([1-prob,prob])*proportion;
            
        case 'uniform'
            proportion = schemeValues{d,2};
            maxN = schemeValues{d,3};
            prob = 1/(maxN+1);
            if(numel(probDist) < maxN)
               probDist(end+1:maxN+1) = 0; 
            end
            probDist = probDist + prob*proportion;
            
        case 'poisson'
            proportion = schemeValues{d,2};
            lambda = schemeValues{d,3};
            useInvCDF = false;
            if(size(schemeValues,2)>=4 && ~isempty(schemeValues{d,4}))
                if(schemeValues{d,4} < 1)
                    useInvCDF = true;
                end
            else
                schemeValues{d,4} = .999;
                useInvCDF = true;
            end
            if(useInvCDF)
                %Uses an inverse CDF of the poisson distribution to
                %calculate the maximum N needed to capture the fraction of
                %the 
                maxN = poissinv(schemeValues{d,4},lambda);
            else 
                maxN = schemeValues{d,4};
            end
            if(numel(probDist) < maxN)
                probDist(end+1:maxN+1) = 0;
            end
            range= 0:maxN;
            prob = lambda.^range.*exp(-lambda)./factorial(range);
            prob = prob/(sum(prob(:)));
            probDist(1:maxN+1) = probDist(1:maxN+1) + prob*proportion;
            
        case 'custom'
            proportion = schemeValues{d,2};
            maxN = numel(schemeValues{d,3});
            prob = schemeValues{d,3}(:)';
            prob = prob/(sum(prob(:)));
            if(numel(probDist) < maxN)
               probDist(end+1:maxN) = 0; 
            end
            probDist = probDist + prob*proportion;
            
    end
end



cumDist = cumsum(probDist(1:end-1));

fluorArray = zeros(numUnits,numel(cumDist));
randVals = rand(numUnits,1);
 
for u = 1:numUnits
    fluorArray(u,:)= randVals(u) >= cumDist;
end

%fluorNumPerUnit
varargout{1} = sum(fluorArray,2);

%fluorIndexIntoUnit
if(nargout>1)
    idxArray = repmat((1:numUnits)' ,1,numel(cumDist));
    fluorArray = fluorArray.*idxArray;
    varargout{2} =  fluorArray(fluorArray ~= 0);
end

%{
for ii=1:numel(cumDist)
    actDist(ii) = numel(varargout{1}(varargout{1} == ii-1));
end
%}

%using hist seems to be an innacurate way to test the created distribution
%[counts,centers]=hist(fluorNumPerUnit,numel(unique(fluorNumPerUnit)));
%[counts2,centers2]=hist(fluorNumPerUnitBinary,numel(unique(fluorNumPerUnitBinary)));

%{
close all;
PrepareFigureSize;
figure;
subPlot(1,1,1);
scatter(0:numel(probDist)-1,probDist);
figure;
subPlot(2,2,1);
scatter(0:numel(actDist)-1,actDist/numUnits);
subPlot(2,2,3);
scatter(0:numel(actDist)-1,actDist/numUnits - probDist(1:numel(actDist)));
subPlot(2,2,2);
scatter(0:numel(actDist2)-1,actDist2/numUnits);
subPlot(2,2,4);
scatter(0:numel(actDist2)-1,actDist2/numUnits - probDist(1:numel(actDist2)));
display(sprintf('The sum method took %d seconds',sumWay));
display(sprintf('The binary method took %d seconds',binaryWay));
1+1;
%}
%figure;
%subPlot(1,1,1);
%scatter(0:numel(actDist)-1,actDist/numUnits,'filled');

