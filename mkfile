all:V: generate solve

clean:V:
	rm -f *.pgm *.ppm sines.txt tmp.log *.png

sines.txt: generate
     ./generate -t 3 -s 7 -n 20

000.pgm: sines.txt
	count=0
	for i in leo*.pgm
    do
		file=$(printf "%03d.pgm\n" $count)
        convert $i -flip -flop $file
        # cp $i $file
        let "count = count + 1"
	done

test:V: 000.pgm solve
    rm -f matches-*.pgm
    time ./solve -t 3 -i 36

graphics/derivative.png: solve
    ./solve -i 70 | 9 tee tmp.log
    cd graphics
    ./derivative.sh

%: %.c helpers.c
	gcc -g $stem.c -fopenmp -lm -o $stem