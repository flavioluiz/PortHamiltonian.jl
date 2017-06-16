using PyCall
using PyPlot
@pyimport matplotlib.animation as anim


#surf(disc.disc_data.flow.xi,t,yout[:,1:N])

plt = PyPlot
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
#ax = plt.axes(xlim=(0, 0.5), ylim=(0.00323439, 0.00323441))
ax = plt.axes(xlim=(0, 0.5), ylim=(0.000, 0.005))
global lines = [ax[:plot]([], [], lw=2)[1];ax[:plot]([], [], lw=2)[1]]

# initialization function: plot the background of each frame
function init()
    global lines
	for line = lines
		line[:set_data]([], [])
	end
	plt.legend(["Nonlinear", "Linear"])
    return (lines,None)
end

# animation function.  This is called sequentially
function animate(i)
	k = 1
	global lines
	for line = lines
		x = [xi]
		if k == 1
			y = [yout[i*4,N+1:end-1][:]]
		else
			y = [youtlin[i*4,N+1:end-1][:]]
		end
		line[:set_data](x, y)	
		#line = plt.plot(x,y)	
		k += 1
	end
    return (lines,None)
end

# call the animator.  blit=True means only re-draw the parts that have changed.
myanim = anim.FuncAnimation(fig, animate, init_func=init,
                               frames=2500, interval=1);

#myanim[:save]("PyPlots-sinplot.mp4", extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])