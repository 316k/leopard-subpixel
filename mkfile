SUBPIXEL=0.5
W=1920
H=1080

.all:V: generate solve subpixel translation subpixel-reference error

reset:V: clean-cam
	rm -f {leo,phase_ref}*.pgm sines.txt tmp.log *.png

clean-cam:V:
    rm -f {phase_cam,matches,0,1,2,3,4,5,6,7,8,9,ref}*.{pgm,ppm}

sines.txt: generate
     ./generate -t 4 -s 3 -n 20 -w $W -h $H

000.pgm: translation sines.txt
    ./add-noise.sh
	# count=0
	# for i in leo*.pgm
    # do
	# 	file=$(printf "%03d.pgm\n" $count)
    #     echo $i to $file
    #     # convert $i -flip -flop $file
    #     # cp $i $file
    #     ./translation -x $SUBPIXEL -y $SUBPIXEL $i > $file
    #     let "count = count + 1"
	# done

ITERATIONS=36

matches-final.ppm: 000.pgm solve
    rm -f matches-*.pgm
    time ./solve -t 4 -i $ITERATIONS
    cp matches-$[$ITERATIONS - 1].ppm matches-final.ppm

ref.ppm: subpixel-reference
    ./subpixel-reference -w $W -h $H -x $SUBPIXEL -y $SUBPIXEL > ref.ppm
    
    for dim in x y
    do
        gnuplot -e "set terminal png size 840,480; set xtics 50; filename='subpixel-$dim-ref'; binwidth=0.001; set xrange [-10:1000]; set output 'histo-ref-$dim.png'" histogram.gnuplot
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
                    ./solve -b $nb_boxes -v $nb_values -t 4 -i 36 > /dev/null
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
        gnuplot -e "set terminal png size 840,480; set xtics 50; filename='subpixel-$dim-vals'; binwidth=0.001; set xrange [-10:1000]; set output 'histo-subpixel-$dim.png'" histogram.gnuplot
    done

graphics/derivative.png: solve
    ./solve -i 70 | 9 tee tmp.log
    cd graphics
    ./derivative.sh

&: &.c helpers.c
	gcc -g $stem.c -fopenmp -lm -o $stem
