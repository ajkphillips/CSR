function s = state_final(t)

global sleeptime

% Define sleep/wake state as a function of time

if sleeptime == -1,

    s = (mod(t,24)>8);
    
elseif sleeptime == 0,
    
    s = (t>8);
    
elseif sleeptime == -2, % Adapted Van Dongen baseline for Wehr fit
    s = (mod(t,24)>7.65);
else
%s = (mod(t,24)>(sleeptime)); % Define s=1 when awake and s=0 when asleep
s = (t<(13*24)).*(1 - (mod(t,24)>(8-sleeptime)).*(mod(t,24)<8)) + (t>=(13*24)).*(mod(t,24)>6.7);
%s = (mod(t,42.85)>(10));

end

end