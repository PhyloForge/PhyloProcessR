# Install PhyloProcessR

PhyloProcessR is an R package. Some functions also call external command-line
programs. Install only the programs that your workflow uses.

The complete environment supports R 4.2 or later. The package itself requires
R 4.0 or later.

## 1. Get the repository

Clone the repository:

```bash
git clone https://github.com/PhyloForge/PhyloProcessR.git
cd PhyloProcessR
```

You can also use **Code > Download ZIP** on GitHub. Extract the ZIP file before
you continue.

## 2. Select an environment method

Use a container when you need the complete reproducible environment. Use Conda
when you need direct access to the installed programs.

### Option A: Use the published container

Use Docker on a local computer:

```bash
docker pull chutter/phyloprocessr:latest
docker run --rm -v "$(pwd)":/app -w /app \
  chutter/phyloprocessr:latest Rscript workflow-5_trimming.R
```

Use Apptainer on a computing cluster:

```bash
apptainer pull phyloprocessr.sif docker://chutter/phyloprocessr:latest
apptainer exec phyloprocessr.sif Rscript workflow-5_trimming.R
```

The Docker Hub account is separate from the GitHub organization. The image
name therefore remains `chutter/phyloprocessr`.

### Option B: Create the Conda environment

Install Miniconda, Miniforge, or another Conda-compatible package manager.
Then create the supplied environment:

```bash
cd setup-files
conda env create -f environment.yml
conda activate PhyloProcessR
cd ..
```

The environment contains the R dependencies and the external programs for the
supplied workflows.

## 3. Install the R package

Install the package from the local checkout when you need an exact local
version:

```r
install.packages("remotes")
remotes::install_local(".", upgrade = "never", dependencies = FALSE)
```

Install the current GitHub version when you do not have a local checkout:

```r
install.packages("remotes")
remotes::install_github(
  "PhyloForge/PhyloProcessR",
  upgrade = "never",
  dependencies = FALSE
)
```

Load the package:

```r
library(PhyloProcessR)
packageVersion("PhyloProcessR")
```

For a reproducible analysis, record the package version or Git commit. Do not
install the moving development branch each time that you run an analysis.

The workflow configuration files set `install.latest.github = FALSE`. Keep
this value for a reproducible analysis. Set it to `TRUE` only when you must test
the newest GitHub code during debugging.

## 4. Check external programs

The workflow scripts check the programs that they use. You can also check a
complete Conda environment in R:

```r
PhyloProcessR::setupCheck(
  anaconda.environment = "/full/path/to/envs/PhyloProcessR"
)
```

If a program is absent, install that program or correct its `*.path` value in
the workflow configuration file.

## 5. Test the installation

Run the small laptop example before you use a new computer or environment:

```bash
cd examples/laptop-example
export PHYLOPROCESSR_BIN="$CONDA_PREFIX/bin"
Rscript run_laptop_example.R
Rscript run_alignment_example.R
Rscript run_trimming_qc_example.R
```

The example uses seven samples and 40 target markers. It does not require a
cluster. See the [example instructions](../../examples/laptop-example/README.md)
for the expected results.

## Troubleshoot Conda write errors

If Conda cannot write package files, add a writable package directory to
`~/.condarc`:

```yaml
pkgs_dirs:
  - /writable/path/conda/pkgs
```

If Conda cannot create environments, add a writable environment directory:

```yaml
envs_dirs:
  - /writable/path/conda/envs
```

Use permanent storage. Conda creates many files.

If dependency resolution fails, create a new environment from the supplied
`environment.yml`. Do not repair a damaged environment by installing packages
one at a time unless you must diagnose a specific conflict.
