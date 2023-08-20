reset
T=10000
set xrange[0:16]
set yrange[0:16]
do for[t=0:T]{
	p 'LJ_sim.txt' i t w p pt 7 ps 1.26 lc 'red'
	set title sprintf("Time t=%d",t)
	#pause 0.01
}
