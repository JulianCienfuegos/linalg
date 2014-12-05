function circle(x,y,r)
% Taken from the mathworks website and modified.
    ang=0:0.01:3*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    for i = [1:length(x)]
        plot(x(i)+xp,y(i)+yp); hold on;
    end
end