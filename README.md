# LipidWrapper 2.0.0

As ever larger and more complex biological systems are modeled in silico, approximating physiological lipid bilayers with simple planar models becomes increasingly unrealistic. When building large-scale models of whole subcellular environments, models of lipid membranes with carefully considered, biologically relevant curvature are essential. We here present a computer program called LipidWrapper, written by Jacob Durrant in the lab of Rommie E. Amaro, capable of creating curved membrane models with geometries derived from various possible sources, both experimental and theoretical. We are hopeful that this utility will be a useful tool for the computational-biology community.

**Note:** This is an unofficial Python 3 update of [the Durrant Lab's LipidWrapper](https://github.com/durrantlab/lipidwrapper). This modernization was largely completed with LLM coding agents, and only superficially validated through light testing, so proceed with caution.

## Installation

```bash
git clone <repository-url> ./LipidWrapper
pip install ./LipidWrapper
```

## Examples

The program download includes an "examples" directory that demonstrate how to use the software. Examples are provided showing how to generate lipid-bilayer models wrapped around equations, PDB point files, DAE models, and image-defined surfaces. Test PDB, DAE, and PNG files are included.

## Command-Line Parameters

### Methods for creating a surface mesh

- `--surface-equation`: Generate a surface mesh from a python-formatted equation defining z, given x and y. The `--min-x`, `--max-x`, `--min-y`, and `--max-y` parameters are used to specify the region over which the function should be evaluated. The `--step-x` and `--step-y` parameters define the x-y distance between adjacent points. Python functions from the math, numpy, and scipy modules can be used. Example: `--surface-equation "z = 250*numpy.sin(x*x/60000 +y*y/60000)"`

- `--surface-filename`: If this parameter specifies a file with the PDB extension, a surface mesh is generated from the coordinates of the PDB atoms. Example: `--surface-filename mymesh.pdb`

- `--surface-filename`: If this parameter specifies a file that does not have the PDB extension, the file is assumed to be a gray-scale image, where black represents regions that are topologically low, and white represents regions that are topologically high. The `--min-x`, `--max-x`, `--min-y`, and `--max-y` parameters are used to specify the region where the mesh should be generated. The `--step-x` and `--step-y` parameters define the x-y distance between adjacent points. The `--max-height` parameter determines the height of the bilayer model at those locations where the image is white; black regions are assigned a height of 0. This feature is only available if the python PIL module has been installed on your system. Example: `--surface-filename mymesh.png`

### The initial lipid model

- `--lipid-pdb-filename`: This parameter specifies a PDB file containing an all-atom model of a planar lipid bilayer. LipidWrapper will wrap this lipid around the user-generated mesh. Example: `--lipid-pdb-filename lipid.pdb`

- `--lipid-headgroup-marker`: A unique atom representing the headgroup of each lipid residue must be specified. The `--lipid-headgroup-marker` accepts a comma-separated lists of atom specifications (RESNAME_ATOMNAME). If either RESNAME or ATOMNAME is omitted, any value will be accepted. By default, LipidWrapper identifies lipid headgroups by looking for any atom named "P" (_P) or any atom named "O3" belonging to a cholesterol molecule (CHL1_O3). Example: `--lipid-headgroup-marker "_P,CHL1_O3"`

### Methods for resolving lipid clashes

- `--delete-clashing-lipids`: It's common for lipids to sterically clash at the interface of two adjacent surface-mesh tessellated triangles. If this flag is set, any clashing lipids are deleted. Example: `--delete-clashing-lipids`

- `--clash-cutoff`: If you do choose to delete clashing lipids, this parameter determines how close two atoms must be (in Angstroms) to constitute a steric clash. Example: `--clash-cutoff 2.0`

- `--fill-holes / --no-fill-holes`: Deleting lipids often leaves holes in the membrane. If enabled, LipidWrapper tries to fill the hole. Default: `--fill-holes`

- `--fill-hole-exhaustiveness`: Essentially, how long LipidWrapper should try to fill the holes. Example: `--fill-hole-exhaustiveness 10`

- `--clashing-potential-margin`: Lipid clashes occur at the edges of adjacent tessellated triangles. If these triangles are very large, it's faster to only check for clashes and holes near the triangle edges. This variable specifies how far from the edges, in Angstroms, that LipidWrapper should look for clashes and holes. Example: `--clashing-potential-margin 25.0`

- `--very-distant-lipids-cutoff`: LipidWrapper determines if two lipids clash by comparing the distance between every atom in the first lipid with every atom in the second lipid. This can be computationally expensive. However, sometimes two lipids are so distant from each other, that it's obvious there are no clashes, making the pair-wise comparison unnecessary. Before performing this expensive pair-wise comparison, LipidWrapper calculates the distance between one atom of each lipid. If this distance is greater than this user-specified cutoff, the program will simply assume there are no clashes. WARNING: Remember to consider the width of your lipid bilayer when choosing this value. Adjacent lipids on opposite sides of the bilayer can seem distant when considering the distance between their headgroups, for example. Example: `--very-distant-lipids-cutoff 50.0`

- `--triangle-center-proximity-cutoff-distance`: Lipid steric clashes/holes typically occur between lipids that belong to adjacent tessellated triangles. However, if tessellated triangles are small enough, clashes are possible between lipids that belong to non-adjacent triangles as well. Consequently, in addition to checking for adjacency, LipidWrapper also checks the distance between the triangle centers, using this user-specified value as a cutoff. Example: `--triangle-center-proximity-cutoff-distance 50.0`

- `--memory-optimization-factor`: When the tessellated triangles are very large and consequently contain many individual lipids, the extensive pairwise distance comparisons required can result in memory errors. This parameter tells lipid Wrapper to divide the list of atoms being compared into smaller chunks. The pairwise distance comparison is performed piecewise on each chunk-chunk pair and so uses less memory, albeit at the expensive of speed. Only increase the value of this parameter if you run into memory errors. Example: `--memory-optimization-factor 1`

### Additional options

- `--number-of-processors`: Using multiple processors can significantly increase the speed of the LipidWrapper algorithm. Example: `--number-of-processors 8`

- `--show-grid-points`: Aside from producing PDB coordinates for lipid atoms, additional coordinates will be appended to the bottom of the output containing "atoms" named "X" that specify the location of the surface mesh points. Example: `--show-grid-points`

- `--create-triangle-tcl-file`: A separate file named "triangles.tcl" will be generated containing a tcl script that can be run in VMD to visualize the mesh surface. Example: `--create-triangle-tcl-file`

- `--output-directory`: If an output directory is specified, all LipidWrapper output files, as well as additional files representing the intermediate steps required to build the final bilayer, will be saved in that directory. Example: `--output-directory ./my_output/`

- `--use-disk-instead-of-memory`: For very large systems, storing the growing model in memory can be problematic. If this flag is set, the growing model will be stored on the hard disk instead. However, expect longer execution times. Example: `--use-disk-instead-of-memory`

- `--compress-output`: Depending on the user options selected, LipidWrapper output can require a lot of disk space. If this flag is set, the output will be automatically compressed using the gzip algorithm (Lempel-Ziv coding LZ77). The files can be uncompressed with the UNIX gunzip utility, or similar Windows-based packages. Example: `--compress-output`

## Example

```bash
lipidwrapper \
  --surface-equation "z = 250*numpy.sin(x*x/60000 +y*y/60000) * (-numpy.sqrt(x*x+y*y)/(560 * numpy.sqrt(2)) + 1)" \
  --min-x 500 --max-x 1000 \
  --min-y 500 --max-y 1000 \
  --step-x 25 --step-y 25 \
  --lipid-pdb-filename lipid.pdb \
  --lipid-headgroup-marker "_P,CHL1_O3" \
  --delete-clashing-lipids \
  --clash-cutoff 1.0 \
  --fill-holes \
  --fill-hole-exhaustiveness 10 \
  > lipid_model.pdb
```
