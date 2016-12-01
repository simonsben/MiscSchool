function [pts] = circ(cent, rad)
    %Function plots circle that can be used for inpolygon function

    x = [linspace(-rad, rad) linspace(rad, -rad)];
    initl = length(x);

    for i = 1:length(x)/2
        y(i) = sqrt(rad^2 - x(i)^2) + cent(2);
    end
    for i = length(x)/2:length(x)
        y(i) = -sqrt(rad^2 - x(i)^2) + cent(2);
    end
    for i = 1:length(x)
        x(i) = x(i) + cent(1);
    end
    
    pts(:, 1) = x';
    pts(:, 2) = y';
end