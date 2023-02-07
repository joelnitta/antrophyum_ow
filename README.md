# antrophyum_ow

Code to run biogeographic analyses for Chen et al. "Systematics and biogeography of the Old World fern genus *Antrophyum*".

All code is in R. Workflow is maintained with the [targets](https://github.com/ropensci/targets) package.

## Dependencies

You must have R and Quarto installed.

This code was developed with R v4.2.1 and Quarto v1.3.14. Other versions may work but are not guaranteed.

## Setup

Place the data files (`Antrophyum.newick`, `Antrophyum_geodata.data`, `manual_dispersal_multipliers.txt`) into `_targets/user/data`.

Start R and run `renv::restore()` to install all needed packages.

## Running the code

Within R, run `targets::tar_make()`. You will see the workflow start running. It may take around 1 hr depending on your computer. When everything is done, a report (`bgb.docx`) and supplemental figures (PDF files) will be output to `reports`.

## License

All code except for `cladistics.csl` is under the [MIT license](LICENSE).

`cladistics.csl` by Steven Vidovic is licensed under a [Creative Commons Attribution-ShareAlike 3.0 License](https://creativecommons.org/licenses/by-sa/3.0/). No changes were made to the file.