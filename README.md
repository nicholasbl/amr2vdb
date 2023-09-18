# amr2vdb
Convert an AMReX mesh to a uniform grid VDB file. This is useful for post-hoc visualization of an AMR-Wind or PeleC run. 

## Installation
This tool uses a buildroot script to download and compile needed packages. This should be as simple as running the `install_deps.py` script.

If the script completes, a `third_party` directory should have been created.

To build `amr2vdb`, you can then 
- `mkdir build`
- `cd build`
- `cmake <path to amr2vdb directory>`
- `make -jN` where N is the number of cores on your machine

When complete, the `amr2vdb` executable should be static, and easily relocatable. 

## Usage

`amr2vdb` takes a TOML file to describe what products should be built. An example is provided in `convert_opts.toml`. Pass your TOML file as the first argument to the tool. Additional arguments are considered to be overrides to the TOML file. So, for example, 

```amr2vdb convert_trim.toml input=\"plt07400\" output=\"amr_struct.vdb\"```

will run the tool, consuming the configuration in the TOML, and will then apply the TOML fragments on the command line to that configuration, in this case, setting the input and output.

Be advised that we are going from an adaptive mesh to a uniform grid. Expect lots of memory usage, and large output files.


## Options

```toml

input = "path" # path to input pltfile
output = "path.vdb" # path to destination vdb

debug = 1 # optional debugging output

amr               # Convert an AMR deck into a vdb, using VDB interpolation
amr.variables     # list of variable names to extract, as strings
amr.blur          # global blur or smoothing options (uses a mean filter)
amr.blur.strategy # one of "none" | "after_every" | "after_last" | "at_end"
                  # after_every => Blur after each interpolation stage as we convert
                  #                a coarse grid to a fine grid
                  # after_last  => Blur after we have collected all coarse grids
                  #                but before we have processed the finest grid
                  # at_end      => Blur final assembled vdb grid
                  
amr.blur.radius   # blur radius
amr.blur.iterations # times to apply blur

amr.save_amr = "name" # save AMR grid fidelity regions as a new grid with name.
                      # useful for understanding the AMR structure

post # postprocessing on VDB grids

post.compute_magvort # if velocity scalars are present, compute magnitude of curl
post.merge_velocity  # combine velocity scalars into a vector valued grid
post.compute_mgd     # compute the magnitude of the gradient of density

post.trim # trim all grids based on the value of a scalar grid

```
