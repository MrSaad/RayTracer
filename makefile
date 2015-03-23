all: 
	gcc -o raytrace rayTracer.c -lm -I/usr/X11R6/include -L/usr/X11R6/lib -lX11 

clean:
	rm raytrace
