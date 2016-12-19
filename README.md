# leopard-subpixel

A leopard-pattern-based computer vision algorithm

## Build

You'll need `gcc` to build these programs

```bash
for i in generate solve subpixel translation subpixel-reference
do
    gcc -g $i.c -fopenmp -lm -o $i
done
```

## Try it

### 1. Run `./generate`

Generates phase shifted leopard patterns.

Use `./generate -help` to see the list of options and default values.


### 2. Capture images

The captured images must be numbered `000.pgm`, `001.pgm`, ...

Alternatively, if `Imagemagick` is installed, you can run `./add-noise.sh` to generate test images


### 3. Run `./solve`

Generates a correspondance map between camera pixels and projected pixels (aka the leopard patterns previously generated).

Use `./solve -help` to see the list of options.


### 4. Run `subpixel`

Use `./subpixel matches-29.ppm` to interpolate and get subpixel accuracy.

This will generate the file `subpixel.ppm`, which should be your final image.
