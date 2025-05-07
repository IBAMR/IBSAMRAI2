# IBSAMRAI2

Backwards-compatible fork of SAMRAI 2.4.4 containing some additional patches,
both for improving performance and for fixing problems with modern compilers.

IBSAMRAI2 assumes that the C++ compiler's default language standard is C++11 or
newer. This is the case for all major compilers released after 2017. If you must
compile IBSAMRAI2 with an older compiler then you must add `-std=c++11` or
equivalent to `CXXFLAGS`. See `./configure --help` for more information in
setting up compilation flags.
