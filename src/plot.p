set terminal gif animate delay 20
set output 'out/heat_2d.gif'
unset border
unset xtics
unset ytics
set cbrange [0.0:0.8]

do for [i=0:99] {
    plot 'out/step_'.(i).".data" with image notitle
}