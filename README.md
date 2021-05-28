# string-amber
Adaptive string method implementation in sander (AmberTools 21.3)

If you use this software in your work, please cite:
https://doi.org/10.1021/acs.jpca.7b10842

### Installation
Copy all the files from sander/ and lapack/ to the corresponding directories in 
`amber20_src/AmberTools/src`. Then set `-DMPI=TRUE` in `amber20_src/build/run_cmake` and 
continue with regular AmberTools installation.

The current files correspond to the version 21.3
of AmberTools. The code *may* work with another subversion, but if you experience
any problems, check the current version by running
```
./update_amber --version
```
Then, if the version is lower than 21.3, run
```
./update_amber --update-to=AmberTools.3
```
and if it's newer, downgrade by running
```
./update_amber --revert-to=AmberTools.3
```

A Dockerfile for AmberTools21+ASM installation in a Docker container can be 
found [here](docker/Dockerfile). AmberTools21.tar.bz2 (as obtained from
https://ambermd.org/GetAmber.php) and [run_cmake](docker/run_cmake) files must 
be available in the build context to build the image.