function h = redondo(X,r,h,forma)
 
    if ~exist('h','var')
      h=figure();
    end
    if ~exist('forma','var')
        figure(h);
        hold on;
        th = 0:pi/50:2*pi;
        if length(r)==1
            xunit = r * cos(th) + X(1);
            yunit = r * sin(th) + X(2);
        else
            xunit = r(1) * cos(th) + X(1);
            yunit = r(2) * sin(th) + X(2);
        end
        plot(xunit, yunit,'LineWidth',2);
        hold off
    else
        figure(h)
        hold on
        th = 0:pi/50:2*pi;
        if length(r)==1
            xunit = r * cos(th) + X(1);
            yunit = r * sin(th) + X(2);
        else
            xunit = r(1) * cos(th) + X(1);
            yunit = r(2) * sin(th) + X(2);
        end
        plot(xunit, yunit ,forma);
        hold off
    end