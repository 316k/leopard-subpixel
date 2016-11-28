.all:V: generate solve enhance subpixel translation

clean:V: clean-cam
	rm -f {leo,phase_ref}*.pgm sines.txt tmp.log *.png

clean-cam:V:
    rm -f {phase_cam,matches,0,1,2,3,4,5,6,7,8,9}*.{pgm,ppm}

sines.txt: generate
     ./generate -t 4 -s 3 -n 20

000.pgm: translation sines.txt
    # ./add-noise.sh
	count=0
	for i in leo*.pgm
    do
		file=$(printf "%03d.pgm\n" $count)
        # convert $i -flip -flop $file
        # cp $i $file
        ./translation -d 0.3 $i > $file
        let "count = count + 1"
	done

matches-final.ppm: 000.pgm solve
    rm -f matches-*.pgm
    time ./solve -t 4 -i 36
    cp matches-35.ppm matches-final.ppm

test:V: matches-final.ppm subpixel
    ./subpixel matches-final.ppm > subpixels
    # Stuff
    echo "bin(x,width)=width*floor(x/width); set terminal png; set output 'subpixel-x.png'; plot 'subpixel-x-vals' using (bin(\$1*100,5)):(1.0) smooth freq with boxes" | gnuplot
    echo "bin(x,width)=width*floor(x/width); set terminal png; set output 'subpixel-y.png'; plot 'subpixel-y-vals' using (bin(\$1*100,5)):(1.0) smooth freq with boxes" | gnuplot
    # error comparaison

graphics/derivative.png: solve
    ./solve -i 70 | 9 tee tmp.log
    cd graphics
    ./derivative.sh

&: &.c helpers.c
	gcc -g $stem.c -fopenmp -lm -o $stem
