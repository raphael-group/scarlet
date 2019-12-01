# SCARLET

SCARLET is an algorithm that reconstructs tumor phylogenies from single-cell DNA sequencing data. SCARLET uses a loss-supported model that constrains mutation losses based on observed copy-number data. 

SCARLET is implemented in Python and uses Gurobi for optimization. 

## Contents
1. [Setup](#setup) 
	- [Dependencies](#dependencies)
	- [Using Gurobi](#usinggurobi)
2. [Running SCARLET](#runningscarlet)
	- [Input](#input)
	- [Output](#output)
	- [Usage](#usage)
3. [Example](#example) 

<a name="setup"></a>
## Setup


<a name="dependencies"></a>
### Dependencies

- Python 2.7 ([anaconda distribution](https://www.anaconda.com/distribution/) recommended)  
- Gurobi 

### SCARLET Setup
To use SCARLET, first clone the SCARLET repository locally. 

```git clone git@github.com:raphael-group/scarlet.git```

SCARLET is then run using Python and doesn't require any compilation. SCARLET does however require an installation and a valid lincense for the Gurobi Optimizier. See the following section for details. 

<a name="usinggurobi"></a>
### Using Gurobi
SCARLET uses the Gurobi optimizer. To setup Gurobi, first download the [Gurobi Optimizer](https://www.gurobi.com/downloads/). 

Every run of SCARLET uses Gurobi, which requires a valid license pointed to by the enviromental variable `GRB_LICENSE_FILE`. Gurobi licenses are freely available for academic users and can be acquired at [https://www.gurobi.com/academia/](https://www.gurobi.com/academia/). There are two types of academic licenses available:

1. **Individual license**. This license can be obtained easily by any academic user with an institutional email. This license is user and machine-specific, meaning that the user needs to acquire a different license for every machine. Assuming the license is stored at `/path/to/gurobi.lic,` the user can easily use it by the following command:

	```export GRB_LICENSE_FILE="/path/to/gurobi.lic"```


2. **Multi-user license**. This license can be used by multiple users on any machine in a cluster. This license can be obtained freely but needs to be requested by the IT staff of the user's institution. This license is typically used in a machine cluster and requires the following command:

	```module load gurobi```

<a name="runningscarlet"></a>
## Running SCARLET

<a name="input"></a>
### Input Files
SCARLET takes as input two files. The first describes the read counts and copy-number profile assignments for each cell. The second describes the copy-number tree and the set of supported losses for each edge in the copy-number tree. 

1. **Read count file**. For each cell, this file describes the copy-number profile assignment and variant and total read counts for *m* mutations. The input file is a comma-separated file, where the first line is the header and rows correspond to cells, with the following format.
	```
	cell_id, c, [mut1]_v, [mut1]_t, [mut2]_v, [mut2]_t, ... [mutm]_v, [mutm]_t
	
	```
	- `cell_id` --- unique identifier for each cell
	- `c` --- [*non-negative integer*] copy-number profile assignment
	- `[muti_v]` --- [*non-negative integer*] the variant read count for mutation *i* in a particular cell
	- `[muti_t]` --- [*non-negative integer*] the integer total read count for mutation *i* in a particular cell. By definition, the total read count must be greater than or equal to the variant read count. 

	Note that in the header `[muti]` may be replaced with unique identifiers for mutations, but the rest of the format must remain as specified.
	

2. **Copy-number tree file**. This file describes the copy-number tree that relates the copy-number profiles of the observed leaves, as well as the set of supported losses for each copy-number transition. This file takes the format of an comma-separated edge list, where the first two columns are integer copy-number profiles, followed by a list of supported losses on that transition. 

	```
	c_1, c_2, [mut1], [mut3], ..., [muti]
	c_2, c_3, [mut2], [mut5], ..., [mutj]
	```
	
	Copy-number profile assignments and mutation names must be consistent with those found in the read count file. 

Examples of these files can be found in the `example/` directory.

<a name="output"></a>
### Output Files

1. **Binary mutation matrix file** (`[prefix].B`). This file describes the binary presence (`1`) or absence (`0`) of each mutation in each cell. This is a comma-separated file where the first line is the header and rows correspond to cells, of the following format. 
	```
	cell_id,[mut1], [mut2], ..., [mutm]
	
	```

2. **Ternary mutation matrix file** (`[prefix].T`). This file is formated in the same way as the binary mutation matrix file but also indicates mutation losses (`2`). That is, if an ancestor of a cell contained a mutation but it was subsequently lost, the mutation matrix entry for that cell will be a `2` instead of a `0`. 

3. **Mutation matrix likelihood file** (`[prefix].LL`). This file contains a single line which gives the log-likelihood of the mutation matrix. 
<a name="usage"></a>
### Usage

SCARLET can be run from the command line as follows.

```
python code/scarlet.py [read count file] [copy-number tree file] [output prefix]
```

<a name="example"></a>
## Example


Example input is provided in the `example` directory. It can be run as 

```
python code/scarlet.py example/read_counts.csv example/tree.csv example/output
```








