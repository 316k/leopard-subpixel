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

### 1. Run `generate`

Generates phase shifted leopard patterns.

Use `./generate -help` to see the list of options.


### 2. Capture images

*TODO* : fill this section.

### 3. Run `solve`

Generates a correspondance map between camera pixels and projected pixels (aka the leopard patterns previously generated).

Use `./solve -help` to see the list of options.


### 4. Run `subpixel`

Use `./subpixel matches-29.ppm` to interpolate and get subpixel accuracy.
This will generate the file `subpixel.ppm`, which should be your final image.
