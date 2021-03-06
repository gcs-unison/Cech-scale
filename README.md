# Cech Scale

This program calculates the cech scale of a system of n disk in 2D. It can also calculate the cech scale for systems of greater dimentions but the system must not have more than 3 disks.

### Prerequisites

You only need a C++ compiler compatible with C++11 standard or above. Also cmake is recommended to easily compile the program.

### Compiling

If you have cmake and are under a GNU/Linux distribution, you can use these commands:

```
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
cd build
make
cd ..
```

And it will compile the program and output the file "Cech-scale" with the program.

## Running

The program by default reads the disk system from the file "texfiles/disks_system.txt"
and writes the results to the file "textfiles/cech_results.txt".

You can also specify others files to read from and write to by calling the program
with the filenames you'd like, e.g for a linux system:

```
./Cech-scale disks.txt results.txt
```

And this would read a disk system from the file "disks.txt" and write the
results of the program to "results.txt"

## Authors

* **Luis Sotomayor** - *Rewriting of the code* - [SanLF](https://github.com/sanlf)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* The code corresponds to the algororithm Cech Scale, which was established in
the paper "A numerical approach for the filtered generalized Cech complex",
by Jesus F. Espinoza, Rosalia Hernandez-Amador, Hector Alfredo Hernandez
Hernandez and Beatriz Ramonetti Valencia.
* Tim Voght's C implementation of the intersection of two circles.

## More info

You can find more info about the original article [here](https://arxiv.org/abs/1809.08175)

A more in depth manual of the project can be found [in the manual directory](https://github.com/gcs-unison/Cech-scale/blob/master/manual/Cech_scale_manual.pdf)

