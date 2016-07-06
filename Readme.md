# TUMmer

TUMmer is a drop-in replacement for MUMmer, being ten times faster on whole genomes. It is based on an enhanced suffix array instead of a suffix array. This makes it much faster, but also requires more memory.

TUMmer does not find all MUMs which overlap in the query. However, these only account for three percent of all MUMs and thus should not affect your analysis.


# Installation

Grab the latest [stable release](https://github.com/kloetzl/TUMmer/releases) and unpack it. TUMmer has [libdivsufsort](https://github.com/y-256/libdivsufsort) as a recommended requirement. It can be installed via the package manager on common Linux systems. If you did get the source, not as a tarball, but straight from the git repository, you will also need the autotools.

Building is as easy as follows.

    $ autoreconf -fi -Im4  # optional when building from tarball
    $ ./configure
    $ make
    $ make install


# Usage

TUMmer follows Unix calling conventions:

    $ tummer foo.fasta
    > AE005674
       1         1        57
      65        65       165
     244       226       166
     411       393       165
     577       559        28

The following options (some with the same functionality as in MUMmer) are supported:

`-b` Compute forward and revere complement matches; default: forward only  
`-j`, `--join` Treat all sequences from one file as a single genome. This might render the position field of the output useless.  
`-l`, `--min-length <INT>` Minimum length of a MUM; uses p-value by default  
`-p <FLOAT>` Significance of a MUM; default: 0.05  
`-r` Compute only reverse complement matches; default: forward only  
`-v`, `--verbose` Prints additional information  
`-h`, `--help` Display help and exit  
`--version` Output version information  


# License

Copyright © 2016 Fabian Klötzl <kloetzl@evolbio.mpg.de>  
License GPLv3+: GNU GPL version 3 or later.

This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. The full license text is available at http://gnu.org/licenses/gpl.html.

Some files may be licensed differently.
