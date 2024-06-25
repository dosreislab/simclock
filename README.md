# Simclock

A simple R package to simulate phylogenies under relaxed clocks.

Prepare a timetree and use the `relaxed.tree` function to generate a phylogeny
with branch lengths in substitutions per site under the desired clock model.
The package implements the geometric Brownian motion (a.k.a correlated rates)
and independent log-normal rates models from Rannala and Yang (2007).

The tree in substitutions per site can then be used with other software
packages (such as Evolver or Seq-Gen) to simulate sequence alignments.

To install the package use
```
devtools::install_github("dosreislab/simclock")
```

## Examples

File `misc/pri10s.tree` contains a timetree for 10 primates species. We will use
this tree to simulate a molecular alignment using the GBM (autocorrelated) rates
model. To do this example, you need to have installed the `ape` R package, and
the Evolver program (from the PAML software package). You also need to know how
to use the command line or terminal. Copy all the files in the `misc/` directory
over to a new directory of your choosing and start R in this new directory.

```
require(ape)
tt <- read.tree("pri10s.tree")

# The primate timetree:
plot(tt)
```

The tree above has branch lengths in millions of years. We will now simulate a
relaxed tree with branch lengths in units of substitutions per site. The mean
mutation rate will be 4e-4 substitutions per site per year, and the diffusion
rate will be 2.6e-3.

```
reltt <- relaxed.tree(tt, model="gbm", r=.04e-2, s2=.26e-2)
plot(reltt)
write.tree(reltt, file="pri10s-relaxed.tree")
```

File `MCbase.dat` has the parameters required by Evolver to simulate the
sequence alignment, except for the tree. Evolver is part of the PAML pakcage for
phylogenetic analysis (Yang 2007). Using your favorite text editor, open this
file. The file will look something like this (see PAML's documentation for full
details on the file's format):

```
 0     * 0,1:seqs or patterns in paml format (mc.paml) format ...
 -1   * random number seed (odd number)

10 1000 1  * <# seqs>  <# nucleotide sites>  <# replicates>
-1         * <tree length, use -1 if tree below has absolute branch lengths>

* SIMULATED RELAXED TREE GOES HERE

7          * model: 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV
0.88892  0.03190  0.00001  0.07102  0.02418 * kappa or rate parameters in model
0.2500  5     * <alpha>  <#categories for discrete gamma>

0.25318  0.32894  0.31196  0.10592    * base frequencies
  T        C        A        G
```

Open file `pri10s-relaxed.tree`, which contains the tree we just simulated, and
copy the tree into the `MCbase.dat` file, so that it looks something like this
(your branch lengths may look different because it is a random simulation):

```
 0     * 0,1:seqs or patterns in paml format (mc.paml) ...
 -1   * random number seed (odd number)

10 1000 1  * <# seqs>  <# nucleotide sites>  <# replicates>
-1         * <tree length, use -1 if tree below has absolute branch lengths>

(Tree_shrew:0.06045474241,((Bushbaby:0.02886311866,Mouse_lemur:0.01663566402):0.007803769466,(Tarsier:0.02621398983,(Marmoset:0.01010700438,(Rhesus:0.01145630454,(Orangutan:0.006589005919,(Gorilla:0.003915870165,(Chimp:0.003059385907,Human:0.003492534071):0.0009415776338):0.003288254269):0.003757859822):0.003802887905):0.009319410718):0.001770574892):0.00930836744);

7          * model: 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV
0.88892  0.03190  0.00001  0.07102  0.02418 * kappa or rate parameters in model
0.2500  5     * <alpha>  <#categories for discrete gamma>

0.25318  0.32894  0.31196  0.10592    * base frequencies
  T        C        A        G
```

In a terminal, while inside the same directory containing the `MCbase.dat` file, type

```
evolver 5 MCbase.dat
```

If the command above does not work, you either don't have PAML installed, or you
don't have the PAML programs in your system's PATH variable (check PAML's
website for details). If successful, Evolver will generate a file called
`mc.paml` containing the simulated nucleotide alignment. In our example above,
there will be 10 simulated alignmnets each with 10 species and of length 1000
nucleotides. Evolver traverses the tree and simulates mutations as the sequences
evolve along the branches.  

The simulated alignments can now be used in downstream analysis. For example,
you can use a program such as `MCMCtree` to infer the divergence times of the 10
species under the appropriate relaxed clock model, and check whether the
inferred times match the true times in the original timetree.  

The package can also simulate trees with correlated rates among loci. This can be
done with the `correlated.trees` function. For example

```
ilnc <- correlated.trees(pri10s, model="iln", r=.04e-2, s2=.1, n=3, corr=0.9)
lapply(ilnc$trees, plot)
```

will simulate three loci with rate correlation of 0.9 (this is a strong
correlation). You can write these trees to a file and use them to simulate three
sequence alignments (one alignment for each tree) but with the proviso that the
aligned loci have evolved in a correlated fashion.

For example `ilnc$trees[[1]]` will print out the first tree to the screen.

```
for (i in 1:3) ape::write.tree(ilnc$trees[[i]], file=paste(i, ".tree", sep=""))
```

This writes the trees to three separate files in Newick format. You can put
these trees in, for example, MCbase.dat to generate the sequence alignments with
Evolver.

## References

* Rannala and Yang (2007) Inferring speciation times under an episodic molecular
clock. Systematic Biology, 56: 453-466.  
* Yang, Z. (2007). PAML 4: Phylogenetic Analysis by Maximum Likelihood. Molecular
Biology and Evolution 24: 1586-1591.  
