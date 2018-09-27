# Cech Scale

This program calculates the cech scale of a system of n disk in 2D. It can also calculate the cech scale for systems of greater dimentions but the system must not have more than 3 disks.

### Prerequisites

You only need a C++ compiler compatible with C++11 standard or above. Also makefile is recommended to easily compile the program.

### Compiling

If you have makefile you only need to run go inside the project folder and run:

```
make
```

## Running

The program by default reads the disk system from the file "texfiles/Disks-system.txt" and writes the results to the file "textfiles/Cech_scale.txt".
You can also specify others files to read from and write to by calling the program with the filenames you'd like, e.g.

```
./cech.exe disks.txt results.txt
```

## Authors

* **Luis Sotomayor** - *Rewriting of the code* - [SanLF](https://github.com/sanlf)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* The code corresponds to the algororithm Cech Scale, which was established in
the paper "A numerical approach for the filtered generalized Cech complex",
by Jesus F. Espinoza, Rosalia Hernandez-Amador, Hector Alfredo Hernandez
Hernandez and Beatriz Ramonetti Valencia.

