function flooredVal = floorS(val, nS)
    pw = ceil(log10(val)); % Find order of magnitude of val.
    res = 10^(pw-nS); % Resolution to round to.
    flooredVal = floor(val/res)*res; % < change floor() to ceil(), for ceiling equivalent.
end