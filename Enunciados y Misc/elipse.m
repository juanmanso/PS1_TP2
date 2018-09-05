function elipse(a,b,x0,y0,color);
	t=-pi:0.01:pi;
	x=x0+a*cos(t);
	y=y0+b*sin(t);
	hold on;
	if color==1
		plot(x,y,'r','linewidth',2)
	elseif color ==2
		plot(x,y,'m','linewidth',2)
	elseif color ==3
		plot(x,y,'k','linewidth',2)
	elseif color ==4
		plot(x,y,'g','linewidth',2)
	end
end
