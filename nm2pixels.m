function pixLocs = nm2pixels(locs,pixSize)

if(~(numel(pixSize) == 1 && isreal(pixSize) && pixSize>0))
    error('pixSize must be a positive real scalar value.');
end

pixLocs = locs/pixSize;

end