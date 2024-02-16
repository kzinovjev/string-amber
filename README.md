> :warning: The current code is incompatible with AmberTools22 or earlier, 
> so make sure you use it with the version 23.  

# string-amber
Adaptive string method implementation in sander (AmberTools 23.05) and pmemd 
(Amber 22.05)

If you use this software in your work, please cite:
https://doi.org/10.1021/acs.jpca.7b10842

### Installation
Copy all the files from sander/ and lapack/ to the corresponding directories in 
`amber22_src/AmberTools/src`. Then set `-DMPI=TRUE` in 
`amber22_src/build/run_cmake` and continue with regular AmberTools installation.

To install ASM for pmemd, copy all the files from pmemd/ directory to 
`amber22_src/src/pmemd/src` and from there run
```
patch -p1 < ASM_amber_22.05.patch
```

The current files correspond to the version 23.05
of AmberTools and 22.05 of Amber. The code *may* work with another subversion, but if you experience
any problems, check the current version by running
```
./update_amber --version
```
Then, if the version is lower than 23.05, run
```
./update_amber --update-to=AmberTools.5
```
and if it's newer, downgrade by running
```
./update_amber --revert-to=AmberTools.5
```

A Dockerfile for AmberTools23+ASM installation in a Docker container can be 
found [here](docker/Dockerfile). AmberTools23.tar.bz2 file (as obtained from
https://ambermd.org/GetAmber.php) must be available in the build context to 
build the image.
