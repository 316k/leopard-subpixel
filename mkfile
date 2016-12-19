SUBPIXEL_X=0.3
SUBPIXEL_Y=0.1
W=1920
H=1080

.all:V: generate solve subpixel translation subpixel-reference error validate dump_pixels

reset:V: clean-cam
	rm -f {leo,phase_ref}*.pgm sines.txt tmp.log *.png

clean-cam:V:
    rm -f {phase_cam,matches,0,1,2,3,4,5,6,7,8,9,ref}*.{pgm,ppm} subpixel-{x,y}-*

sines.txt: generate
     ./generate -s 3 -n 20 -w $W -h $H

000.pgm: translation sines.txt
    # ./add-noise.sh
	count=0
	for i in leo*.pgm
    do
		file=$(printf "%03d.pgm\n" $count)
        echo $i to $file
        # convert $i -flip -flop $file
        # cp $i $file
        ./translation -x $SUBPIXEL_X -y $SUBPIXEL_Y $i > $file
        let "count = count + 1"
	done

ITERATIONS=66

matches-final.ppm: 000.pgm solve
    rm -f matches-*.pgm
    time ./solve -i $ITERATIONS
    cp matches-$[$ITERATIONS - 1].ppm matches-final.ppm

ref.ppm: subpixel-reference
    ./subpixel-reference -w $W -h $H -x $SUBPIXEL_X -y $SUBPIXEL_Y > ref.ppm
    
    for dim in x y
    do
        sort subpixel-$dim-ref | awk '{ printf("%.3f\n", $1) }' | uniq -c > h
        cat h
        echo no
        gnuplot -e "set terminal png size 840,480; set xtics 50; filename='h'; set output 'histo-ref-$dim.png'" histogram.gnuplot
        
        # gnuplot -e "set terminal png size 840,480; set xtics 50; filename='subpixel-$dim-ref'; binwidth=0.001; set xrange [-10:1000]; set output 'histo-ref-$dim.png'" histogram.gnuplot
    done
       
benchmark-solve: 000.pgm solve error
    for iteration in $(seq 10)
    do
        echo == Iteration $iteration ==
        for nb_boxes in $(seq 7 16)
        do
            for nb_values in $(seq 7 16)
            do
                if test $(echo $nb_boxes ^ $nb_values | bc) -lt 134217728
                then
                    ./solve -b $nb_boxes -v $nb_values -i 36 > /dev/null
                    for i in 5 4 3 2 1
                    do
                        echo -n "$i $nb_boxes $nb_values $(echo $nb_boxes ^ $nb_values | bc) "
                        ./error -l $i.1 matches-35.ppm ref.ppm | grep bad | wc -l
                    done
                fi
            done
        done
    done

test:V: matches-final.ppm subpixel ref.ppm error
    ./subpixel matches-final.ppm
    
    # Cute graphics
    
    for dim in x y
    do
        sort -k 2 subpixel-$dim-vals | awk '{ printf("%.3f\n", $2) }' | uniq -c | sort -k 2 > h
        gnuplot -e "set terminal png size 840,480; set xtics 50; filename='h'; set output 'histo-subpixel-$dim.png'" histogram.gnuplot
    done
    
    echo -n ''>h-proportions
    for i in $(seq 0 3)
    do
        echo -n "$i "
        grep "^$i" subpixel-x-vals | wc -l
    done | awk '{ print $2, " ", $1 }'>h-proportions
    
    gnuplot -e "set terminal png size 840,480; set xtics 50; filename='h-proportions'; set output 'histo-proportions.png'" histogram.gnuplot

presentation.pdf: presentation.md
	pandoc -t beamer -V theme:Berkeley -V colortheme:dolphin presentation.md -o presentation.pdf

&: &.c helpers.c
	gcc -g $stem.c -fopenmp -lm -o $stem
