Authors: Michele Busby (code and specification)

What does it do?
Process a bam file to get all of the alignments for a random set of reads at the read level, not alignment.

What would I want to do that?
There are many ways to randomly downsample a bam. Use this one if:
-You need to set the number of reads beforehand (it won't be exact - it will be some slop, like around a poisson distribution)
-You have multiply mapped reads in your sample
The random downsampling in Picard gives you a random set of *alignments*.  Reads that map to multiple places may or may not have all of their alignments sampled, and will be surely overrepresented in the final sampling. This gets a set of random set by reads, regardless of their alignment.  If in doubt you probably want this one,not the Picard one.

Don't use this one if:
You have huge bam files and nothing is multiply mapped. In that case, you could fill up your memory and it will crash. In that case, use samtools or picard or whatever you like.

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

How to get it to compile:

Compiling the file may or may not be a big hassle because you need bamtools. My make file may NOT be doing the linking quite right. If you know how to fix that please let me know! But the program should still run.

First, though, copy all the files in here to a single directory on your unix server. The server should have g++ is installed (it probably is).


Bamtools

To get this to compile you will need to download Derek Barnett's bamtools API and install it in a folder somewhere. It is here: https://github.com/pezmaster31/bamtools You need the API, not the command line program though is quite useful to have. Then, after you install it, 

--> you need to edit the Makefile in this folder so that everywhere it says "FolderWhereBamToolsIs" you put in the folder where bamtools is located.

Compiling

Go to the folder where everything is installed in type "make". Ignore the warnings about Handy doing stupid stuff. If nothing says ERROR and it makes an executable called ComplexityByStartPos you should be all set.

License:
You can use it.
