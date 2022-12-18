> :warning: The current code is incompatible with AmberTools21, so make sure you use
> it with the version 22.  

> :warning: A string interpolation bug was introduced in May 2022, but the fix wasn't 
> released to the main branch. If you get nonsensical values in the last node in 
> `0_final.string`/`0_final_CV.string`, please, reinstall ASM using the current 
> version of the code.

# string-amber
Adaptive string method implementation in sander (AmberTools 22.3)

If you use this software in your work, please cite:
https://doi.org/10.1021/acs.jpca.7b10842

### Installation
Copy all the files from sander/ and lapack/ to the corresponding directories in 
`amber20_src/AmberTools/src`. Then set `-DMPI=TRUE` in `amber20_src/build/run_cmake` and 
continue with regular AmberTools installation.

The current files correspond to the version 22.3
of AmberTools. The code *may* work with another subversion, but if you experience
any problems, check the current version by running
```
./update_amber --version
```
Then, if the version is lower than 22.3, run
```
./update_amber --update-to=AmberTools.3
```
and if it's newer, downgrade by running
```
./update_amber --revert-to=AmberTools.3
```

A Dockerfile for AmberTools22+ASM installation in a Docker container can be 
found [here](docker/Dockerfile). AmberTools22.tar.bz2 file (as obtained from
https://ambermd.org/GetAmber.php) must be available in the build context to 
build the image.