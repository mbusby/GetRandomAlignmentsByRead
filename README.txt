Authors: Michele Busby (code and specification)

What does it do?
Process a bam file to get all of the alignments for a random set of reads.

What would I want to do that?
The random downsampling in Picard gives you a random set of *alignments*.  Reads that map to multiple places
may or may not have all of their alignments sampled, and will be surely overrepresented in the final sampling.
This gets a set of random set by reads, regardless of their alignment.  If in doubt you probably want this one,
not the Picard one.

How does it work? 
Type:
./GetRandomByRead -bam /seq/picard/flowcell/lane/Solexa-IC/name.bam -out output.bam -nReads 1000000 -seed 49

Required parameters:
-bam Name of the bam file
-out Name of the output file
-nReads The number of reads you want.
Optional:
-seed A random number seed. This is so the program reports the same sampleing every time, or you can give it different
seeds if you want two different random datasets.  Default is my son's birthday.

This uses a lot of memory so you may run out.
If that happens your computer will say "segmentation fault".  You can try:
use LSF
ish
To get to a node with more memory or you may have to submit it to the LSF queue.  See Michele for more help.


It crashed. I am sad. And I can't really figure out how it works.:
The utilities in this directory are quick and dirty so they aren't as full tested and not 
optimized for usability and don't have proper documentation.  (It's not you, it's me.)  So if you're 
stuck just contact Michele at mbusby@broadinstitute.org or stop by.