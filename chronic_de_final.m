function [T,Y] = chronic_de_final(xi,params,T)

[T,Y] = rk4(@(t,y) odeChronic(t,y,params),0,T,xi,T*10);

function dy = odeChronic(t,y,params)

    % dy = change in variables
    % y(1) = Atot = concentration of adenosine molecules
    % y(2) = Rtot = concentration of adenosine receptors

    Ab = 0.5*((y(1)+y(2)+params{1}/(1-params{5})) - sqrt((y(1)+y(2)+params{1}/(1-params{5}))^2 - 4*y(1)*y(2)));
    Au = y(1) - Ab - params{5}; % concentration of unbound adenosine molecules 
    Ru = y(2) - Ab; % concentration of unbound adenosine receptors
    
    dy = [state_final(t)*(1/params{8})*(params{7}-y(1)) + (1-state_final(t))*(1/params{9})*(params{6}-y(1));
        params{10}*(Ab - y(2)*params{4})];
        

end % end ode

end