# hirshfeld

[![Build Status](https://travis-ci.org/lidaobing/hirshfeld.png?branch=master)](https://travis-ci.org/lidaobing/hirshfeld)

hirshfeld can help you calculate the hirshfeld charge from the
gaussian's fchk file. hirshfeld charge is defined in
doi:10.1007/BF00549096 .

following is an example for calculation hirshfeld charge.

```
$ hirshfeld HCN.fchk
No.     Atomic  elctron         charge
1       1         0.871666      +0.128334
2       6         5.947401      +0.052599
3       7         7.180643      -0.180643
```

I write this program because gaussian03's hirshfeld charge module is
buggy and not free. This program is licensed under GPL.

## DOWNLOAD

you can download the stable version from
http://code.google.com/p/hirshfeld/ .


## INSTALL

run

```
$ ./cofigure; make
$ sudo make install
```

## USAGE

run

```
$ hirshfeld foo.fchk
```

foo.fchk is convert from gaussian's chk file by

```
$ formchk foo.chk foo.fchk
```

## CUSTOMIZE

hirshfeld need the electro density information for each atom, I only
provide H, C, N, O, P, S's data. these data is obtained by gaussian
using ub3lyp/6-311++G(d,p) and the multiplicity with lowest energy. If
you need atoms whose data does not included in this package or you don't
like the data provided me, you can obtain data file by following steps:
(for example, you need the data of P)

1. do the calculation of a single atom and keep the chk file (named
15.chk), pay attention to select the correct multiplicity.

2. run

```
$ formchk 15.chk 15.fchk
$ mkdir -p ~/.local/share/hirshfeld
$ /usr/local/lib/hirshfeld/convert 15.chk \
  > ~/.local/share/hirshfeld/15.data
```

The directory of the data file is followed by XDG Base Directory
Specification[1](http://freedesktop.org/wiki/Specifications_2fbasedir_2dspec) (current version is 0.6). By default, you can put the
data file in one of these directories:

* $HOME/.local/share/hirshfeld/
* /usr/local/share/hirshfeld/
* /usr/share/hirshfeld/

it also influeced by the environment variable XDG_DATA_DIRS, if you need
more information, check [1](http://freedesktop.org/wiki/Specifications_2fbasedir_2dspec).

## DEVELOPMENT VERSION

you can obtain the development version by
`git clone git://github.com/lidaobing/hirshfeld.git`

## BUGS

fill bug report or wishlist at
https://github.com/lidaobing/hirshfeld/issues

## THANKS

First thanks to F. L. Hirshfeld, who introduce Hirshfeld charge to the
world, maybe you want to cite the paper by F. L. Hirshfeld:

Hirshfeld, F.
Bonded-atom fragments for describing molecular charge densities
Theoret. Chim. Acta (Berl.), 1977, 44, 129-138

Also thanks to Prof. WU Yundong, I wrote this program when I were in his
group.
