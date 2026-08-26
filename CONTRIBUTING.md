# Contributing to PhyloProcessR

Thank you for helping improve PhyloProcessR. Bug reports, documentation fixes,
tests, and focused code contributions are welcome.

## Before opening an issue

- Search the existing issues to avoid duplicates.
- For a bug, include the PhyloProcessR version or Git commit, operating system,
  R version, relevant workflow and configuration, and the smallest reproducible
  example you can provide.
- Include the names and versions of external executables when the problem occurs
  in a workflow stage that invokes them.
- Remove credentials, private sample identifiers, and unpublished sequence data
  from logs and examples.

Questions about scientific defaults or changes to analytical behavior should be
opened as a discussion issue before implementation. Describe the biological
rationale, expected result, and any reference implementation or publication.

## Development setup

Fork and clone the repository, then install the package dependencies. The full
container or Conda environment is needed only when exercising workflows that
invoke external bioinformatics programs.

```R
install.packages(c("devtools", "testthat"))
devtools::install_dev_deps()
devtools::load_all()
```

Run the fast unit-test suite with:

```R
devtools::test()
```

Before submitting a pull request, regenerate documentation and run a package
check:

```R
devtools::document()
devtools::check(args = c("--no-manual", "--no-build-vignettes"))
```

## Pull requests

- Keep each pull request focused on one change.
- Add or update tests for changed behavior.
- Update roxygen comments and regenerate `man/` and `NAMESPACE` when the public
  API changes. Do not hand-edit generated `.Rd` files.
- Preserve existing scientific behavior unless the change fixes a demonstrated
  bug or has been discussed with the maintainer.
- Do not include generated analysis output, large sequence files, credentials,
  or machine-specific paths.
- Describe any external-tool versions and small fixture data used for workflow
  integration testing.

By participating, you agree to follow the project
[`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md).
