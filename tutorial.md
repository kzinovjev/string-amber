# Adaptive string method tutorial

Table of contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Preparation](#preparation)
  - [CVs file](#cvs-file)
  - [guess file](#guess-file)
  - [STRING file](#string-file)
- [Running the job](#running-the-job)
- [String preparation stage](#string-preparation-stage)
- [String optimization stage](#string-optimization-stage)
  - [convergence.dat](#convergencedat)
  - [STOP_STRING](#stop_string)
  - [*.string](#string)
  - [\#.dat](#dat)
  - [*.PMF](#pmf)
  - [*_CV.PMF](#_cvpmf)
  - [force_constants.dat and node_positions.dat](#force_constantsdat-and-node_positionsdat)
  - [*.REX](#rex)
- [Umbrella Sampling stage](#umbrella-sampling-stage)
  - [0_final.string and 0_final_CV.string](#0_finalstring-and-0_final_cvstring)
  - [final_parameters.dat](#final_parametersdat)
  - [*_final.PMF](#_finalpmf)
  - [*_CV_final.PMF](#_cv_finalpmf)  
  - [\#_final.dat](#_finaldat)
- [Analysis](#analysis)
  - [reorder_trj.py](#reorder_trjpy)
  - [get_pmf_cv_values.py](#get_pmf_cv_valuespy)
  - [get_ts_frames.py](#get_ts_framespy)
- [Other features](#other-features)
  - [Keep the system close to the path during the US stage](#keep-the-system-close-to-the-path-during-the-us-stage)
  - [Restart the PMF calculation for a converged string](#restart-the-pmf-calculation-for-a-converged-string)
- [Final remarks](#final-remarks)

Introduction
------------

This tutorial explains how to set up, run, and analyze an ASM calculation for the reaction between chloromethane and chlorine ion:

IMAGE

This is one of the examples from the [original ASM paper](https://pubs.acs.org/doi/10.1021/acs.jpca.7b10842). The system is small enough for the entire protocol to finish in about 15 minutes on a single cluster node (20-30 CPUs). Alternatively, a complete setup with all the input and output files can be found [here](). 

The string will be optimized in the space of 4 collective variables (CVs): two carbon-chlorine (C-Cl) distances, the hybridization of the carbon atom (measured as the signed point-plane distance between the carbon and the plane formed by the hydrogens), and Cl-C-Cl angle. The latter is not really needed to trace the reaction process (and was not included in the calculation presented in the paper) but is included here to illustrate that the method works well when mixing the CVs of very different nature (distances and angles).

Prerequisites
-------------
To run the preparation and analysis scripts make sure you have access to Python3 installation with libraries `molmod`, `mdtraj` and `parmed` available. All can be installed using `pip`.
 
Preparation
-----------

To run any ASM calculation one needs the [amber topology of the system]() (inpcrd/parm7 file) and the [reactants]() and [products]() structures equilibrated using the same potential that will be used during the ASM calculation. It is also important for the potential to be exactly the same during the relaxation of reactants and products. For instance, in this tutorial, semiparabolic restraints on C-Cl distances were used to keep the system from dissociating. The [restraint definition]() is identical in both relaxations and will be kept the same for the ASM calculation. 
 
It is also beneficial (but not mandatory) to have short unbiased MD trajectories of [reactants]() and [products]() - these can be used to obtain the positions of the minimum free energy path (MFEP) endpoints and thus accelerate the string convergence. This approach is used here, as explained in the next section.
 
 
Apart from the system topology and the structures, ASM calculation requires multiple input files that define all the aspects of the calculation (definitions of the CVs, initial guess etc.). The file format is not very user friendly, but for the most error-prone parts, there are scripts available that make the file generation easier.
 
#### *CVs* file
The `CVs` file, as the name suggests, contains the definitions of the CVs used. The easiest way to generate it is to use the [generate_cvs.py](utils/generate_cvs.py) script:
```bash
$ generate_cvs.py <topology> cvs.def > CVs
```
where `cvs.def` is a simplified CVs definitions file: 
```
BOND :1@C :1@CL1
BOND :1@C :2Cl-
PPLANE :1@C :1@H2 :1@H1 :1@H3
ANGLE :1@CL1 :1@C :2@Cl-
```
The format of `cvs.def` is much simpler compared to the resulting `CVs` file. Each line describes a single CV and contains the keyword defining the type of the CV, followed by ambermasks for the atoms involved in the CV. Fours types are supported:

1. `BOND` - interatomic distance (2 atoms)
2. `ANGLE` - bond angle (3 atoms)
3. `DIHEDRAL` - dihedral angle (4 atoms)
4. `PPLANE` - signed point plane distance (4 atoms)

For `PPLANE`, the first atom defines the point, while the remaining three atoms define the plane. Note that it is a signed distance, so the sign of the resulting CV depends on the order of atoms defining the plane.

#### *guess* file
The `guess` file contains the initial guess for the MFEP. it is defined as a set of points in the CV space with the following format:
```
<number of guess points>
<CV_1 in point 1> <CV_2 in point 1> ... <CV_D in point 1>
<CV_1 in point 2> <CV_2 in point 2> ... <CV_D in point 2>
...
<CV_1 in point N> <CV_2 in point N> ... <CV_D in point N>  
```
Here D and N are the number of CVs and the number of guess points, respectively. The first and the last points are especially important since they define the ends of the MFEP (the reactants and products free energy minima) in the chosen CV space. While it is possible to put approximate values for these points and let the algorithm optimize them, the default behavior is to keep the ends fixed during the string optimization. This dramatically improves the convergence compared to a "free" string. The positions of the minima have to be obtained by taking the average values of the CVs along unbiased MD trajectories of reactants and products. This can be done with [trj_to_cvs.py]() script:
 ```bash
$ trj_to_cvs.py CVs <trajectory> <topology> <first frame> <last frame> 
 ```
It calculates the values of the CVs (defined in the `CVs` file) for each frame of the provided trajectory within the specified range. If the first and the last frames are not specified, the entire trajectory is used. If only the first frame is specified, the frames until the end of the trajectory will be parsed. After printing the values for individual frames, the script also prints the average values of the CVs over the specified range. These values are the ones that should be used in the `guess` file. 
 
So, the entire `guess` file can be constructed from the trajectories provided [here]() as follows:
 ```bash
$ echo 2 > guess
$ trj_to_cvs.py CVs react.nc parm7 41 | tail -n 1 >> guess
$ trj_to_cvs.py CVs prod.nc parm7 41 | tail -n 1 >> guess
 ```
Here we only take the last 10 frames out of 50 to calculate the averages. The following `guess` file is generated:
```
2
  1.76510e+00  2.99074e+00 -3.48126e-01  1.50092e+02
  2.97581e+00  1.76765e+00  3.13631e-01  1.50445e+02
```
Note that although the system is symmetric, there are slight differences between the values in the reactants and products due to statistical error. However, the differences are small enough to be ignored.

#### *STRING* file
The last file that has to be generated manually is the `STRING` file. It is the main configuration file and its presence instructs sander to actually perform the ASM calculation. There are lots of parameters that can be defined in the `STRING` file, however, in the vast majority of cases the default values work well, so nothing has to be specified. In this example, all the parameters will be left at their default values, the only two things that have to be provided are the output directory for the results of the calculation and the path to the `guess` file:
```
$STRINGSETUP
dir = "results/"
guess_file = "guess"
$END
```
It always makes sense to put all the ASM output files in a separate directory, because a large number of files is generated and keeping them together with the sander outputs (which are also quite a few) makes it even messier.

Running the job
---------------
The ASM code relies on multisander for job parallelization - instead of running a single MD simulation, a set of independent simulations are performed, with sander command-line arguments for each one being specified in the so-called "groupfile". This allows specifying different initial structures and output files for each simulation. Here the first half of the nodes will be initiated from reactants and the rest will be initiated from the products. Also, it is useful to provide separate sander input files for each node, with random number seed (*ig* parameter) explicitly specified. In case of any problems, this allows rerunning exactly the same simulation, but, for example, writing out the trajectory frames more often. Both input files and the groupfile can be generated using the following script ([in.sh]()):   
```bash
#!/bin/bash

# If groupfile exists, remove it
if [ -f string.groupfile ];
then
  rm string.groupfile
fi

PARM=../system/ch3cl.parm7  # path to topology file
REACT=../system/react.rst   # path to reactants structure
PROD=../system/prod.rst     # path to products structure
SEED=12345                  # Some random number 

# Number of string nodes is provided as command line argument
nodes=$1
for i in `seq 1 $nodes`; do
  # generate sander input for node $i 
  sed "s/__SEED__/$((SEED + i))/g" in > $i.in
  
  # use reactants structure for the first half of the nodes and products for the rest
  if [ $i -le $(($nodes/2)) ];
  then
    crd=$REACT
  else
    crd=$PROD
  fi
  
  # write out the sander command for node $i to the groupfile
  echo "-O -rem 0 -i $i.in -o $i.out -c $crd -r $i.ncrst -x $i.nc -inf $i.mdinfo -p $PARM" >> string.groupfile
done
```
Note that in principle multiple different reactants and products structures can be used to initialize the nodes, but here only one of each is used for simplicity. Also, note that the `*.in` files are generated based on the template file `in`. The only difference between the `in` file and a regular sander input is the following line:
```
ig=__SEED__
```
The `__SEED__` string is a placeholder that is replaced by a different number for each sander input. The total number of steps (`nstlim`) should be set to a significantly large number for the job not to die before the string converged and enough sampling was acquired for PMFs. The idea is that the job should always be killed manually by the user, after checking that no more sampling is needed. So, it is better to put some unreasonably large number here.

With both the `groupfile` and all the `*.in` files in place, sander can be executed as:
```bash
$ mpiexec sander.MPI -ng 24 -groupfile string.groupfile
```
Here -ng specify the number of groups (string nodes). While in the example below 24 string nodes are used, it should also work well with any number between 20 and 30. The following command might need to be executed before running sander to ensure compatibility of the code with newer versions of MPI libraries:
```bash
export I_MPI_COMPATIBILITY=3
```
Apart from that, ASM calculations are in no way different from any other MPI job. For convenience, example scripts for [PBS]() and [SLURM]() are also provided.

String preparation stage
------------------------
Immediately after the job starts and all the input files are parsed, the ASM calculation enters the "preparation" stage. By default, it lasts 1 ps and during this stage the force constant of the string biasing potential (that keeps the simulation close to the position of the string node in the CV space) gets gradually increased from 0 to the target value. So, in each string node, the system is slowly brought to the target region of conformational space, without causing numeric instabilities. If nevertheless, some nodes become unstable (which often happens in the vicinity of the TS) the number of steps used for preparation can be increased using `preparation_steps` parameter in the `STRING` file.

String optimization stage
-------------------------
After the preparation is done, the string optimization stage starts. That when the string starts to move towards the MFEP and both the force constants and node positions along the path get optimized. Multiple files are generated during this stage and described below.

#### *convergence.dat*

This is arguably the most important file - it tracks the convergence of the string optimization by measuring the RMSD between the current state of the string and the entire history of the string states. The first column of the file is the timestep and the second is the RMSD value. When the RMSD plateaus it indicates that the calculation reached convergence and the string now fluctuates around the MFEP:

```
gnuplot> plot 'convergence.dat' w l
```
IMAGE


#### *STOP_STRING*
The `STOP_STRING` file is not created by the code, but should be created by the user to indicate that the string had converged at it is time for the code to switch to the Umbrella Sampling stage. The file contains only two lines: the optimization step at which the user assumes that the string had already converged and another step, that defines the window used to take averages of the string node positions in the CV space and all the parameters needed to define the pathCV and run Umbrella Sampling simulations. For example, from the plot above it seems that the string had already converged after 10000 steps of optimization. However, we want to wait until step 15000 to calculate the necessary averages over the last 5000 steps. Then, `STOP_STRING` will look as follows:
```
10000
15000
```
Once the `STOP_STRING` is created, the code will recognize it, wait until step 15000, take the averages and switch to the Umbrella Sampling stage.

In principle, just knowing how to interpret `convergence.dat` and create `STOP_STRING` is sufficient to use the method as a "black box". As described in the next section, once the simulation enters the Umbrella Sampling stage, there is nothing the user has to do - the Umbrella Sampling along the pathCV will be performed and the PMFs will be calculated fully automatically. However, when something goes wrong (and it certainly will) the rest of the files generated during the preparation stage provide lots of valuable information for debugging, so it is helpful to know how to interpret them. 

#### **.string*
The `*.string` files contain the state of the string at a given step of the optimization. By default these are written every 100 steps, it can be changed with `output_period` parameter in `STRING`. The format is straightforward: each column corresponds to a CV and each row to a string node. For example, the file `0.string` (which is written even before the string optimization starts) contains the initial guess, which is just a linear interpolation between the two points provided in the `guess` file:
```
gnuplot> plot for [i=1:4] '0.string' u i w l title ''.i
``` 

IMAGE

#### *\#.dat*
These files contain the values of the CVs in each node during the string optimization. So, the number in column 3 of line 100 in file `10.dat` is the value of the third CV (the carbon hybridization) in string node 10 at step number 100. Note that the steps are counted from the beginning of the preparation stage, not from the beginning of the whole simulation. The remaining columns contain technical information that can be ignored.

#### *\*.PMF*
Contain approximate PMFs calculated, by default, using the last 2 ps of sampling and are printed every 2 ps. The frequency can be changed by specifying the desired number of steps with `buffer_size` parameter in the `STRING` file (the sampling from all the nodes is stored in a buffer and after writing out the PMF the buffer is emptied, hence the name of the parameter). Note that these PMFs are **VERY**  approximate and should only be used to get a qualitative idea of whether the calculation behaves as expected. The contents of all the columns are described in the file header. For instance, the values of the RC are given in column 1, while the PMF is found in column 6:
```
gnuplot> plot '10000.PMF' u 1:6 w l
```

IMAGE

#### **_CV.PMF*
Contain the decomposition of the PMF in CV components according to Eq. 33 in [ASM paper](https://pubs.acs.org/doi/10.1021/acs.jpca.7b10842). Since the decomposition introduces additional uncertainty, this information is even noisier than the PMFs and is only kept for consistency with the *_final_CV.PMF files that are written out during Umbrella Sampling (see below).

#### *force_constants.dat* and *node_positions.dat*
These files contain the time evolution of the force constants and positions of the string nodes along the MFEP (rescaled to [0:1] interval), respectively. While not as important as `convergence.dat`, this information can also be used to judge, whether the simulation converged or to diagnose any convergence problems. For instance, weird behavior of node positions (rapid changes, several nodes located in the same place) indicates that too few string nodes are used, while wild jumps in the force constants might be a symptom of some structural instability. Since the example used here is very simple, both the force constants and node positions change smoothly:

```
gnuplot> plot for [i=1:24] 'force_constants.dat' u i w l title ''
```
IMAGE

```
gnuplot> plot for [i=1:24] 'node_positions.dat' u i w l title ''
```

IMAGE

It takes some practice to tell a "good" plot from a problematic one, so for beginner users just looking at `convergence.dat` to decide when to stop the string optimization should be enough.

#### **.REX*
All the `*.REX` files contain information about the replica exchange procedure:
- `exchanges.REX`: list of all exchanges. The first column is the step, the second - index of the first of the nodes involved in the exchange. So, if there is a line in `exchanges.REX` containing `2000   12`, then at step 2000 of string optimization, there was an exchange between nodes 12 and 13. This file is rarely used.
- `histogram.REX`: each line of this file contains the total number of exchanges in the corresponding node. For example, line 10 contains the number of exchanges between nodes 10 and 11. It is preferable to have at least several exchanges for each of the nodes. If there are zeroes in this file, it might mean that too few nodes were used.
- `plot.REX`: allows for a visual representation of the exchanges. The first column is the step and the rest are the indices of the string nodes sampled by each MD simulation:

```
gnuplot> plot for [i=2:25] 'plot.REX' u 1:i w l lw 2 title ''
```

IMAGE

The plot shows how each trajectory "walks" through different string nodes.

Umbrella Sampling stage
-----------------------
Once the number of string optimization steps reaches the second value provided in `STOP_STRING` file (see above), the code calculates MFEP by averaging and interpolating the positions of the string nodes. This provides the definition of the pathCV, which is used as the reaction coordinate for US calculations. The force constants and reference values for US windows are obtained by averaging the data in `force_constants.dat` and `node_positions.dat`. From this point, no user action is required, apart from checking the PMF files (see below) to decide whether enough sampling is acquired and the job can be terminated.

All the files generated during Umbrella Sampling contain "final" in the filename. For instance, all the replica exchange related files have the same meaning as the ones described above, but are called `*_final.REX`. Other files have some differences compared to their counterparts from string optimization and are described below. 


#### *0_final.string* and *0_final_CV.string*
These are the only `*.string` files written out during the US. Both contain the averaged string that was used to define the path CV. The `0_final.string` also contains some technical information used to [restart the US calculation](#restart-the-pmf-calculation-for-a-converged-string). The contents of `0_final_CV.string` are also slightly different: the first column gives the value of the path CV, that corresponds to each US window, while the rest of the columns are values of the CVs, as in the regular `*.string` files described above. This file is therefore useful to relate the PMF with the regions in the CV space. For example, the CV values from the window sampling the maximum of the PMF approximately describe the TS geometry.

#### *final_parameters.dat*
Contains the positions and force constants of the US windows. Used to [restart US calculation]().

#### **_final.PMF*
The most important files - contain the actual PMF along the pathCV. In contrast to the PMF files generated during the string optimization, these are integrated using the entire sampling obtained during the US stage. So, the file `10000_final.PMF` contains the PMF calculated based on 10000 steps of US. The more sampling is acquired, the smaller is the statistical error, which is given in the column 7. The rule of thumb is to stop finish the calculation when the error at the TS is < 1 kcal/mol.

#### **_CV_final.PMF*
Same as for the string optimization, these files contain the CV-decomposition of the PMF. The interpretation of the values can be found in the [ASM paper]().

#### *\#_final.dat*
As the `#.dat` files from string optimization, the `#_final.dat` files contain values of the CVs in the corresponding window, but they also contain the values of the path CVs - *s* (the progress along the path, **the** reaction coordinate) and *z* (distance from the path) as well as the parameters of the US window. The format is the following:
```
<force constant along s> <position along s> <force constant along z> 0.00000E+00
<s in step 1> <z in step 1> <CV_1 in step 1> ... <CV_D in step 1> ... (technical stuff)
...
```
The first line contains the force constants and RC values that define the US bias. Note that the bias along *z* coordinate is always centered at 0, because this coordinate is used to keep the system close to the path. Also, note that by default the force constant along *z* is also 0 - the system is free to move in the directions orthogonal to the path. The following lines contain the values of the path CVs and the values of the CVs defined in the `CVs` file as well as some technical stuff not important in most of the cases.

Analysis
--------
The most important structural features impacting the reaction (C-Cl distances and CH3 hybridization) are included in the CV list, they are available in the `#_final.dat` files together with the path CV values and can be analysed directly. However, if we, for example, would like to measure the average chlorine-chlorine distance at the TS, it is not so easy: in contrast to a simple RC, like the antisymmetric transfer coordinate, the path CV can not be trivially calculated from a structure in VMD or PyMOL. Also, if we just want to look at how the structures in the vicinity of the TS look like, it can not be easily done: because of the replica exchange, the MD trajectories written out by sander are not in one-to-one correspondence with the string nodes and US windows. To solve these issues, several postprocessing scripts are available and described below. All the scripts can be found [here]().


#### reorder_trj.py
Takes the MD trajectories from ASM calculation and reorders the frames so that the new trajectories do correspond to the string nodes during the string optimization and US windows during the US stage. It must be executed from the working directory (the one containing the `STRING` file) with a single argument - the system topology: 
```bash
$ reorder_trj.py parm7
```
The script automatically locates the results directory and used the information from `plot.REX` and `plot_final.REX` to rearrange the frames into new trajectories. The results are written into `reorder_trj_results` directory. Separate trajectories are generated for string optimization (`string` subdirectory) and US (`pmf` subdirectory) stages.

#### get_pmf_cv_values.py
Once the trajectories are reordered, it is possible to extract the values of the path CV for each frame from the `#_final.dat` files. This is done by running [get_pmf_cv_values.py]() script from the working directory without any arguments:
```bash
$ get_pmf_cv_values.py
```
It will create a set of `#.dat` files in `reorder_trj_results/pmf` with path CV values for each frame of the reordered trajectories. The format of `#.dat` files is exactly as of `#_final.dat`, so it includes *s* CV, *z* CV, and the CVs used for string optimization.

#### get_ts_frames.py
A common task is to extract and analyze only the frames corresponding to the TS of the reaction. While it can be done by using the converged PMF and the results from [get_pmf_cv_values.py](), it can be done much faster with [get_ts_frames.py]() script. It also runs from the working directory and requires two arguments: the topology and the PMF file which will be used to determine the position of the TS. For example:
```bash
$ get_ts_frames.py parm7 --pmf results/20000_final.PMF
```
The script will generate a new trajectory named `ts.nc` by default (can be changed with `--output` argument) that contains uniquely the frames belonging to the TS. By default, all the frames with path CV within 0.05 of the PMF maximum are extracted. This threshold can be changed with `--thr` argument. 

Other features
--------------
#### Keep the system close to the path during the US stage
By default, when switching from the string optimization to the US, the bias that restrains the movement of the system in the directions orthogonal to the path is removed. It is done because it is assumed that when the string reaches the MFEP, the system will not be able to drift away from it because of the underlying free energy landscape. However, it is often the case that the MFEP is not well defined, so the simulation still can explore regions far from it. It might cause numerical instabilities and artifacts in the PMF. In this case, the orthogonal bias should be kept in place by changing the `remove_z_bias` parameter (set to `.true.` by default in the `STRING` file:
```
remove_z_bias = .false.
```
By doing so, the harmonic bias over *z* path CV (the one that measures the distance from the path) is added during the US optimization.

#### Restart the PMF calculation for a converged string
When the string takes unexpectedly long to converge, it might happen that there is not much time left to acquire enough sampling during the US stage (e.g. due to job time limits on the cluster). In this case, it is possible to take advantage of the converged string and start a new jb directly with US stage, skipping the preparation and string optimization. To do so:
1. Add `only_PMF = .true.`  to the `STRING` file
2. Copy the files `final_parameters.dat` and `0_final.string` from the old results directory to the new **working** directory (the one with the `STRING` file).
3. Get initial structures for all the nodes from the previous job. Keep in mind that the frames in the sander MD trajectories are not properly ordered, so you first should run `reorder_trj.py` script to get them in the correct order. Then, for example, you can simply take the last frames of the `reorder_trj_results/pmf/#.nc` trajectories.
4. Change `in.sh` script to provide the initial structures obtained in step 3 to corresponding windows. For example, if the initial structures are stored in the `init` subfolder of the working directory, the line in `in.sh` that writes to the groupfile can be modified as follows:
```
...  
echo "-O -rem 0 -i $i.in -o $i.out -c init/$i.ncrst $crd -r $i.ncrst -x $i.nc -inf $i.mdinfo -p $PARM" >> string.groupfile
...
```
When running the job prepared that way, no `STOP_STRING` file is needed - it will immediately start generating all the `*final*` files.

Final remarks
-------------
While the method is made as "black box" as possible, in the sense that there are very few parameters to tweak and input files are fairly easy to prepare, there are several things to keep in mind.

First, it seems that the string optimization and especially the path CV are more sensitive to any instabilities compared to regular US simulations. If you struggle with protons flying around or weird bonds being formed during a regular QM/MM US simulation, ASM will probably make it worse. So, it is quite common to add multiple restraints to reinforce the bonds that are not suppose to break etc. 

Second, while ASM can treat very complex multiple-step reactions in a single job, it is a good idea to split a "long" string into "shorter" ones whenever it is possible. So, if you know that your reaction has a stable intermediate, it makes sense to run two separate jobs: *reactants &rarr; intermediate* and *intermediate &rarr; products*. 

Finally, the method works best when there is little friction along the MFEP due to the environment and the rate of the process is mostly limited by a high potential energy barrier that has to be surpassed going from reactants to products. This is exactly the case for chemical reactions, but not so much for more "diffusive" processes, such as large conformational changes or dissociation events. While there is some evidence that ASM can perform well even for diffusive processes, it requires tweaking the parameters controlling the string optimization, which goes out of the scope of this tutorial.