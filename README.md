# FNTexScan

This is companion code for the paper [Mapping and Probing Froggatt-Nielsen Solutions to the Quark Flavor Puzzle](https://arxiv.org/abs/2306.08026), arXiv:2306.08026.

This code is written in the [Nim](https://nim-lang.org) programming language, with some Python code (Jupyter notebooks) for visualization.

Dependencies are stored in the FNTexScan.nimble file.  To install all dependencies, run `nimble check` in the top-level directory.  Then to build the main binary, use:
```
nim c -d:danger FNTexScan/FNTexScan
```

Note that the `nim.cfg` file in that directory is set for building on MacOS with Apple silicon; you may need to adjust it for your machine.

Aside from some additional Nim libraries, the [GNU Scientific Library](https://www.gnu.org/software/gsl/) is a required dependency.  The code must be linked to the GSL using the library flag `-lgslcblas`; see the `nim.cfg` file.  On my machine, I also had to add the library directory to my `$DYLD_LIBRARY_PATH`; since I use homebrew, this amounted to adding the flag

`export DYLD_LIBRARY_PATH=/opt/homebrew/lib:$DYLD_LIBRARY_PATH`

to my `~/.zshrc` file.