# Frequency of Different Amino Acid Sequences in Protein-Protein Interactions

## Project Overview
This project involved creating an algorithm that counted the amount of protein-protein interactions where a sequence of $k$ neighbouring amino acids were involved within a cut-off distance of $\theta$\AA.

This algorithm was deployed on a large data set of 870 PDB files to test the hypothesis that some sequences of amino acids are more frequently involved in protein-protein interactions.

This algorithm was capable of being scaled up for use in a distributed memory parallelised computing environment with hyper-threading. This was done in order to deliver results in a shorter amount of time while working on large data sets.

After experimentally determining the parallel computing configuration that led to the optimal speed-up and efficiency, a final analysis was conducted for all combinations of $k=6$ and $\theta=\{4,5,6\}$\AA.

Using the optimal configuration of six nodes, each with two threads, on the data set of 870 PDB files, led to a speed-up of three times compered to serial computing. A parallel efficiency of 40\% was achieved and the results were computed in a mean time of 13.61 hours with a standard deviation of 3.22 hours.

After analysing the results, it appeared that some sequences of amino acids like ASN-ASN-TYR-ALA-ASP-PHE and ASN-TYR-ALA-ASP-PHE-ASP were indeed involved more frequently in protein-protein interactions. However, even the sequences that were involved in most interactions were never involved in more than 0.5\% and 1.9\% of the total amount of interactions for $\theta=6$\AA and $\theta=10$\AA, respectively. 

Therefore, although some amino acid sequences might appear more frequently in protein-protein interactions, they cannot be used as reliable predictors of sequences that will be involved in most interactions.

## Backgroud
Proteins are macromolecules that influence cell behaviour in many living organisms. These proteins are made up of smaller molecules called amino acids. These amino acids are in turn made up of atoms.

The amino acids that make up a protein molecule are chained together in a linear sequence. This linear sequence of amino acids is referred to as the primary structure of a protein.

However, proteins fold into complex three-dimensional shapes due to the chemical and physical properties of their constituent amino acids and do not resemble their primary structures. This three-dimensional structure of a folded protein is called its tertiary structure.

Protein-protein interactions largely determine the functions of cells and are even used in drug-development. For the purpose of this project, locations of protein-protein interactions are defined as locations (in the tertiary structure) where two sections of one or more proteins are in close proximity.

More specifically, these interactions are defined as locations where the centroid of one amino acid comes within a distance of $\theta$ \r{a}ngstr\"om (\AA) of the centroid of another amino acid. These two amino acids can be on different protein chains or, if they are on the same chain, it will only be considered an interaction if they are separated by at least 20 amino acids in the primary structure.

Centroids of amino acids are defined as the average location of their constituent atoms which are assumed to have point locations. The centroid $(\bar{x}, \bar{y}, \bar{z})$ of an amino acid can be calculated using equation \ref{eq:centroid}.
\begin{equation}
\label{eq:centroid}
    (\bar{x}, \bar{y}, \bar{z}) = \left( \frac{1}{n}\sum_{i=0}^{n-1} x_i, \frac{1}{n}\sum_{i=0}^{n-1} y_i, \frac{1}{n}\sum_{i=0}^{n-1} z_i \right)
\end{equation}
Where $x_i$, $y_i$, and $z_i$ are the $x$, $y$, and $z$ coordinates of the $i\textsuperscript{th}$ atom respectively.

The aim of this project was to test the hypothesis that certain sequences of $k$ consecutive amino acids are more likely to be involved in protein-protein interactions than others.

## File Structure

|---code
|    |
|    |---run_test.py  // the python code that finds interacting amino acid sequences
|    |---test_script.sc  // an axample of the script that was submitted to SLURM to run experiments/jobs

|---singularity_image
|    |
|    |---my_proj_def.def // singularity definition file that will create an environment where my results can be reproduced
|                        // first create singlarity container by running: [sudo] singularity build project.sif my_proj_def.def
|                        // then test by running: singularity exec project.sif mpiexec -np <nb_nodes> python run_test.py <directory/containing/pdb/files> <k> <theta>
|
|---report.pdf  // the project report including deailed methodulogy and reslts
