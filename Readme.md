# TUMmer

TUMmer is a drop-in replacement for MUMmer, being ten times faster on whole genomes.

TUMmer does not find MUMs which overlap in the query. However, these only account for three percent of all MUMs and thus should not affect your analysis.


# Installation

Grab the latest [stable release](https://github.com/kloetzl/TUMmer/releases) and unpack it. TUMmer has [libdivsufsort](https://github.com/y-256/libdivsufsort) as a recommended requirement. It can be installed via the package manager on common Linux systems. If you did get the source, not as a tarball, but straight from the git repository, you will also need the autotools.

Building is as easy as follows.

    $ autoreconf -fi -Im4  # optional when building from tarball
    $ ./configure
    $ make
    $ make install


# License

Copyright © 2016 Fabian Klötzl <kloetzl@evolbio.mpg.de>  
License GPLv3+: GNU GPL version 3 or later.

This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. The full license text is available at http://gnu.org/licenses/gpl.html.

Some files may be licensed differently.
