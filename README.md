# python_modules

Python helper modules and functions are collected globally here.

The intended use is to try to do any computation with a function and put the script in this module.
Then we can use those modules around any analysis we have without copying and pasting code around.
To distribute the analysis we only need to ship this module along with the analysis scripts, install it with pip as shown above and we are ready to reproduce it.


## Installation

First clone this repo from git.
As it includes a data_mc_plotter fork as a submodule, we have to clone recursively, as described [here](http://stackoverflow.com/questions/3796927/how-to-git-clone-including-submodules):

```bash
git clone --recursive git@github.com:jwerthebach/python_modules.git
```

If the repo was already normally cloned (non recursively) then we can add the submodule later with:

```bash
git clone git@github.com:jwerthebach/python_modules.git
# < Forgot to clone recursively... >
cd jopymods/data_mc_plotter
git submodule update --init --recursive
```

In the main repo we can list all submodules by simply using

```bash
git submodule
```

which will list all submodules and their paths.

The repo is a python module folder and can be installed with pip:

```bash
cd python_modules
pip install --user -e .
```

which puts it in the pip site-packages folder and makes it available systemwide if the site-packages are in the path
.
The `-e` option (developer mode) only links the package and let's you change code without reinstalling it over and over again.
This solution avoid manually adding paths to the PYTHONPATH environment variable.

### Uninstall

The module can be uninstalled normally by using `pip uninstall anapymods`.


```bash
pip install --user -e ./jopymods
```

This will make the folder available systemwide and puts it in the pip site-packages folder.
The `-e` option (developer mode) only links the package and let's you change code withour reinstalling it over and over again.
This solution avois manually adding paths to the PYTHONPATH environment variable.
The folder name in `python_modules` defines also the module name under which it gets installed, here `jopymods` (JOhannesPYthonMODuleS).

The intention is to put every function needed more than one time in this module. This avoids a lot of copying and pasting.
To distribute the analysis we only need to ship this module along with the analysis scripts, install it with pip as shown above and we are ready to reproduce it.

Submodules get named after their purposes, e.g.

	plotting
	data_mc_plotter
	...

Put every new category in a seperate folder and add it to the `setup.py` under `packages`.
For some more detail on `setup.py` files see [example here](https://github.com/pypa/sampleproject/blob/master/setup.py).

The module can also be uninstalled normally by using `pip uninstall jopymods`.


## Infos

- In `setup.py` add packages manually or use `find_package()`:

  ```python
  # packages = find_package(where=".")
  packages = [
  	"jopymods",
  	"jopymods.submodule1",
  	"jopymods.submodule2",
  	# and so on
  	]
  ```
- Each subfolder needs an `__init__.py`. Here we can import submodules to create a usability like in e.g. `scipy`. The third line is used to hide dfilenames and imports from the namespace.

  ```python
  from ._tests import *
  from ._sampling import *
  __all__ = [_s for _s in dir() if not _s.startswith('_')]
  ```
- Then we can import stuff with `import jopymods.statistics as amps` and then the `amps` object holds alle the functions.
- Further differentiation is done via subfiles. E.g. the submodule (folder) `statistics` contains `_tests.py` for statistical tests and `sampling.py` for sampling functions. The filenames are prepended with an underscore to hide it from the namespace. The same goes for imports inside the files. Using `import numpy as _np` hides `_np` from the namespace.

## attgen

All functions in the module `attgen` needs to be run in an `icecube` environment. There for the *IceCube Offline Software* must be installed and the environment must be loaded by running the `env-shell.sh`.

## data_mc_plotter

data_mc_plotter is forked from https://github.com/mbrner/data_mc_plotter.git in https://www.github.com/jwerthebach/data_mc_plotter and added as a submodule to this repo.
What follows are short instructions how to deal with this combination.

### fork and setup as submodule

See [github](https://help.github.com/articles/fork-a-repo/) on how to fork a repo and [github](https://help.github.com/articles/syncing-a-fork/) to sync this repo.

Then as described [on git-scm](https://git-scm.com/book/de/v1/Git-Tools-Submodule) we make a submodule from the fork with 

```bash 
git submodule add git@github.com:jwerthebach/data_mc_plotter.git jopymods/data_mc_plotter
```

Then we can add and commit outside like normal.

For changes made inside the data_mc_plotter submodule we need to be inside it.
Let's add a `__init__.py` to make data_mc_plotter availabe within our `jopymods` module:

```bash
cd jopymods/data_mc_plotter
touch __init__.py
git add __init__.py
git commit -m "Made data_mc_plotter work inside jopymods by adding __init__.py"
git push
```

This will push the changes to the data_mc_plotter fork and does nothing with the outside repo.

Outside we can track the change from data_mc_plotter by running `git status` normally.
We will see, that data_mc_plotter is noted as changed `modified:   jopymods/data_mc_plotter (new commits)`.
So we track the data_mc_plotter commits from outside with our normal workflow:

```bash
git add jopymods/data_mc_plotter
git commit -m "Made data_mc_plotter work inside jopymods by adding __init__.py"
git push
```

When cloning the repo somewhere else, the newest submodule revision needs additionally be checked out by calling `git submodule update` from `./python_modules`.

### setup fork

To be able to pull newest changes from the original data_mc_plotter we need to set the upstream to mbrner github with 

```bash
git remote add upstream https://github.com/mbrner/data_mc_plotter.git
```

Thgen syncing works as described in [github help](https://help.github.com/articles/syncing-a-fork):

```bash
git fetch upstream
git checkout master
git merge upstream/master
```

The structure of this repo and most of the `README.md` were taken from [mennthor](https://github.com/mennthor/python_modules).