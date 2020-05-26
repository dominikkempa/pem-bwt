pEM-BWT - Parallel EM construction of BWT from SA
=================================================


Description
-----------

pEM-BWT is an implementation of the parallel external-memory
algorithm to construct Burrows-Wheeler transform from the text
and its suffix array. The algorithm was described in the paper.

    @inproceedings{kk16esa,
      author =    {Juha K{\"{a}}rkk{\"{a}}inen and Dominik Kempa},
      title =     {Faster External Memory {LCP} Array Construction},
      booktitle = {24th Annual European Symposium on Algorithms (ESA 2016)},
      pages     = {61:1--61:16},
      year      = {2016},
      doi       = {10.4230/LIPIcs.ESA.2016.61},
    }

The latest version of pEM-BWT is available from
https://github.com/dkempa/pem-bwt.



Requirements
------------

pEM-BWT has no external dependencies (libraries, cmake, etc).
It only requires:
- g++ compiler supporting the -std=c++0x flag (all modern versions)
- A 64-bit operating system. The current version has been tested
  on Linux/PC.



Compilation and usage
---------------------

The package contains a single Makefile in the main directory. Type
`make` to build the two executables that allow computing the BWT
of a given file sequentially and in parallel. For usage instructions,
run the programs without any arguments. pEM-BWT relies on the prior
computation of the suffix array for the input text. The currently
fastest external memory algorithms for suffix array construction are
pSAscan and fSAIS. Both algorithms are available from
https://github.com/dkempa/.

### Example

The simplest usage for pEM-BWT is as follows. Suppose the text
is located in `/data/input.txt`, the suffix array of `input.txt` is
located in `/data/input.txt.sa5` and is encoded using 40-bit integers.
Then, to compute the BWT of `input.txt` using pEM-BWT, type:


    $ ./compute_bwt_parallel /data/input.txt


This will write the output BWT to `/data/input.txt.bwt`. By default,
the algorithm uses 3.5GiB of RAM and assumes that the input text is
over byte alphabet. A more advanced usage is demonstrated below.


    $ ./compute_bwt_parallel ./input.txt -s ~/data/input.txt.sa -o ../input.txt.bwt -c 4 -i 6 -m 8gi


Explanation:
- The -s flag allows specifying the location and filename of the
  suffix array. The default location and filename is the same as
  input text, with the appended ".saX" suffix, where X is the used
  integer size (see the explanation of -i flag below).
- The -o flag allows specifying the location and filename of the
  output BWT. The default location and filename is the same as input
  text, with the appended ".bwt" suffix.
- The -c flag allows specifying the size of text symbol (in bytes).
  The default size of text symbol is 1 byte. In this example, the
  symbol size is set to 4 bytes. Currently supported are values from
  the range [1, 8].
- The -i flag allows specifying the integer size (in bytes) used to
  encode the input suffix array. The default integer size is 5 bytes.
  In this example, the type is set to a 6-byte integer. Currently
  supported are values from the range [4, 8].
- The -m flag allows specifying the amount of RAM used during the
  computation (in bytes). In this example, the RAM limit is set to
  8gi = 8 * 2^30 bytes (see the explanation below).

Notes:
- The symbol in BWT at position corresponding to suffix of the
  text starting at position 0 (i.e., position ISA[0]) is 0.
- The argument of the -m flag (RAM used during the computation)
  can be specified either explicitly or using common suffixes
  such as K, M, G, T, Ki, Mi, Gi, Ti, which respectively correspond
  to multipliers: 10^3, 10^6, 10^9, 10^12, 2^10, 2^20, 2^30, 2^40.
  Suffix names are not case-sensitive, e.g., Ti = ti, k = K.
- The flags specifying integer type, output filename, etc. can be
  given in any order.
- Filenames passed as command-line arguments can be given as relative
  paths, e.g., `../input.txt.bwt` and `~/data/input.txt.sa` are valid
  paths, see also example above.
- To enable additional statistics about the computation (alternative
  counter of I/O volume and tracing of the disk usage), uncomment line
  with AUX_DISK_FLAGS in the Makefile. When this flag is enabled, the
  computation could slow down thus this flag is disabled by default.
- The usage of the sequential version of the algorithm is identical to
  the parallel version. The only difference in the execution is that
  the sequential version uses exactly one thread for computation.



Troubleshooting
---------------

1. I am getting an error about the exceeded number of opened files.

Solution: The error is caused by the operating system imposing a
limit on the maximum number of files opened by a program. The limit
can be increased with the `ulimit -n newlimit` command. However, in
Linux the limit cannot be increased beyond the so-called "hard limit",
which is usually only few times larger. Furthermore, this is a
temporary solution that needs to repeated every time a new session is
started. To increase the limits permanently, edit (as a root) the file
`/etc/security/limits.conf` and add the following lines at the end
(including the asterisks):


    * soft nofile 128000
    * hard nofile 128000


This increases the limit to 128000 (use larger values if necessary).
The new limits apply (check with `ulimit -n`) after starting new session.



Limitations
-----------

- At present the only limitation in the usage of the algorithm is the
  need to ensure that the limit for the number of opened files in
  the system is sufficiently large to prevent the above error. This
  technical shortcoming will be eliminated in the future versions of
  pEM-BWT.



Terms of use
------------

pEM-BWT is released under the MIT/X11 license. See the file LICENCE
for more details. If you use this code, please cite the paper mentioned
above.



Authors
-------

The authors of pEM-BWT are:
- [Dominik Kempa](https://scholar.google.com/citations?user=r0Kn9IUAAAAJ)
- [Juha Karkkainen](https://scholar.google.com/citations?user=oZepo1cAAAAJ)
