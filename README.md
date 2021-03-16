# string-amber
Adaptive string method implementation in sander (AmberTools19+)

If you use this software in your work, please cite:
https://doi.org/10.1021/acs.jpca.7b10842

### Installation
1. Compile serial version of AmberTools:
   ```bash
   $ ./configure gnu # For GNU compilers
   $ make install # Use -j option to parallelize
   ```

2. Compile parallel version of AmberTools:
   ```bash
   $ ./configure -mpi gnu # For GNU compilers
   $ make install # Use -j option to parallelize
   ```

3. Introduce the changes from the repo to AmberTools/src. Since your AmberTools subversion is most probably different, it is strongly recommended not to copypaste files `sander/Makefile`, `lapack/Makefile`, `sander/qm2_extern_module.F90`, `sander/sander.F90` and `sander/force.F90`, but compare the contents with the files you have and manually introduce only the changes relevant to the string method. Other files can be copied as is.

4. Update the depend file:
   ```bash
   $ ./makedepend > depend # From AmberTools/src/sander directory
   ```

5. Recompile AmberTools without reconfiguring:
   ```bash
   $ make # From AmberTools/src directory
   ```
