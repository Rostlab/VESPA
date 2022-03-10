<details>
  <summary markdown="span"> <u><b>Table of Contents</u></b> </summary>

[[_TOC_]]

</details>

# VESPA -- **V**ariant **E**ffect **S**core **P**rediction without **A**lignments

**VESPA** is a simple, yet powerful Single Amino Acid Variant (SAV) effect predictor based on embeddings of the Protein Language Model [ProtT5](https://github.com/agemagician/ProtTrans).

The single-sequence-based SAV effect prediction is set up in a multistage pipeline that includes (1) generating ProtT5 embeddings, (2) extracting per-residue conservation predictions, (3) (optionally) extracting per-variant log odds ratios, and (4) predicting the effect of all possible amino acid substitutions. Step (4) can be completed by either using **VESPA** with (2) and (3) as input, or by using the computationally more efficient method VESPA-light (**VESPAl**) with only step (2) as input for a small drop in prediction performance.

The specifics of **VESPA** and **VESPAl** can be found in our paper, [_Embeddings from protein language models predict conservation and variant effects_ (Marquet et al., 2021)](https://doi.org/10.1007/s00439-021-02411-y). The performance of **VESPA** when evaluated against SOTA methods can be seen below.
<div align="center"><img src="images/point_main_0823.png" width="75%" height="75%"/></div>

## Precomputed VESPA and VESPAl predictions

Precomputed **VESPA** and **VESPAl** predictions are currently available for 39 _DMS experiments_ here:  [_Supplementary file 3 of (Marquet et al., 2021_)](https://doi.org/10.1007/s00439-021-02411-y).

Furthermore, **VESPAl** predictions for the _human proteome_ (downloaded 22/01/17) with the corresponding FASTA file are available here: [_https://doi.org/10.5281/zenodo.5905863_](https://doi.org/10.5281/zenodo.5905863).

## Usage

The preferred method to install **VESPA** is via pip:

```bash
pip install vespa
```

### Input Files

**Required**: A single `fasta` file containing all your wildtype sequences (note: this file can contain any number of sequences).

{+Optional+}: If you are only interested in a **subset of possible mutations (specific mutations)**, you can add `-m mutations.txt` to the code line in [Quickstart](#Quickstart) (note: per default all mutations are considered). Click [`here`](#mutations-file) to head to the file format explanation of `mutations.txt` directly .

-----

For simplicity of this guide, we will assume a folder containing all data: f.e., the `fasta` file is placed at `data/sequences.fasta` and the (optional) `mutations.txt` at `data/mutations.txt`.

### Quickstart

After installing the repository you can run the following:

```bash
vespa data/sequences.fasta --prott5_weights_cache data/t5_xl_weights
```

This will generate a new folder `vespa_run_directory` in your current working directory. Within this folder you will find three files containing input features (embeddings, conservation prediction, log-odds ratios) and an `output` folder with .csv files containing **VESPAl** predictions for each sequence with all possible mutations.

**WARNING** Creating embeddings requires a powerful GPU! For more details, see [Step 1 on generating embeddings](#step-1-extracting-prott5-embeddings) and [Step 3 on extracting log-odds ratios](#step-3-log-odds-ratio-of-masked-marginal-probabilities).

{+Optional+}: If you want to have the best possible predictions and have a powerful GPU available you can generate **VESPA** predictions by adding `--vespa` to the code line above (per default only **VESPAl** is run). Please note that running **VESPA** will require substantially more time and resources than running **VESPAl**.

Below you can find information on how to run each step of **VESPA** and **VESPAl** individually. We also introduce several optional arguments such as predicting specific mutations instead of the entire mutational landscape.

### The substeps of VESPA and VESPAl

The following steps can be run individually and for any number of sequences contained in a FASTA file. Running the `vespa` script will automatically run all of the below for you. Using the substeps is only required if you are interested in any particular intermediate result.

#### Step 1: Extracting ProtT5 embeddings

To run **VESPA** and obtain SAV predictions, you will need the `protT5` embeddings of your sequences. These you can generate using the [embed.protein.properties](https://embed.protein.properties/) server (make sure to select `ProtTrans T5 XL U50 (ProtT5)`).

Alternatively, if you have a powerful GPU (we recommend at least 12GB of VRAM), you can use the included embedding script on your own machine to generate your protein embeddings.

```bash
vespa_emb --prott5_weights_cache data/t5_xl_weights -o data/embeddings.h5 data/sequences.fasta
```

#### Step 2: Conservation Prediction

**VESPA** and **VESPAl** take per-residue conservation predictions as input. To generate them, run the following (in the `VESPA` folder):

```bash
vespa_conspred data/embeddings.h5 -o data/conspred
```

This will generate 9-state conservation probabilities (per-residue) needed as input for the models. For more details on the conservation prediction please see the paper mentioned above.

{+Optional+}: In case you are interested in generating a file that contains only the predicted conservation classes instead of the assigned class probabilities, add a`--output_classes` flag.

#### Step 3: Log odds ratio of masked marginal probabilities

**WARNING** This step requires a powerful GPU; we recommend at least 12GB of VRAM!

We provide two versions of **VESPA**: **VESPA** and **VESPAl**. They differ in their
predictive performance but also in the required input. If you are only interested to run
**VESPAl**, you can skip this step.

To generate the log odds ratios of masked marginal probabilities required for **VESPA**, run (in the `VESPA` folder):

```bash
vespa_logodds -o data/logodds.h5 data/sequences.fasta
```

**WARNING** This will compute the log odds ratio for every possible SAV of a sequence and thus can take a long time. In case you are only interested in specific mutations, especially when generating output for a large number of sequences or long sequences, we recommend to include a _mutations file_. The file format is described below.

{+Optional+}: To run the log odds script including the mutations file, add `-m data/mutations.txt` to the code line above.

**WARNING** Per default the command above generates the [`.h5` file](#h5-files) required for the subsequent steps.

{+ Optional+}: We provide two options to output a human readable version of the log odds scores by adding `--single_csv data/single_logodds.csv` or `--csv_dir data/csv_dir/` at the end of the code line above. The first option generates a single csv file with all sequences and SAVs. For large sets, we recommend to use the second option, which outputs multiple csv files separated by sequence ID of the given FASTA file.
The format of the csv files is described [below](#log-odds-ratio-output).

#### Step 4: Run **VESPA** and/or **VESPAl**

Now you have all the data required to run **VESPA** and/or **VESPAl**. To explicitly disable one model, add `--no-vespa/ --no-vespal`. To generate predictions, execute the following code (in the `VESPA` folder):

- **Both** (if you computed the conservation prediction and the logodds):

    ```bash
    vespa_run --vespa --vespal -v data/conspred_probs.h5 data/sequences.fasta --T5_input data/logodds.h5 -m data/mutations.txt --output predictions/
    ```

- **Only VESPA** (if you computed the conservation prediction and the logodds):

    ```bash
    vespa_run -v --no-vespal data/conspred_probs.h5 data/sequences.fasta -T5_input data/logodds.h5 -m data/mutations.txt --output predictions/
    ```

- **Only VESPAl** (if you only computed the conservation prediction):

    ```bash
    vespa_run -v --no-vespa data/conspred_probs.h5 data/sequences.fasta -m data/mutations.txt --output predictions/
    ```

**Note:** Running VESPA automatically generates predictions for both models.
The format of the results file is described below at [VESPA and VESPAl output](#vespa-and-vespal-output).

### Additional Information

#### Extracting raw reconstruction probabilities

You might be interested in extracting the raw reconstruction probabilities for each mutation position from T5. These raw reconstruction probabilities are used to calculate the log odds ratio and is explained in more detail in the corresponding publication.
To do so, use:

```bash
poetry run vespa_logodds -r --reconstruction_output data/reconstruction_probas -o data/logodds.h5 -m data/mutations.txt data/sequences.fasta
```

The generated datasets in the [`.h5` file](#h5-files) will contain -1 if a particular mutation was not computed. Specifically positions that were not considered (i.e. the position was not present in the optional `mutation.txt`) will be represented by a vector containing only -1. Otherwise the data contains contain probability vectors that determine the reconstruction probabilities for all amino acids sorted according to the `MUTATION_ORDER` in `config.py`.

#### File Specifications

This section describes a few relevant file formats we use for **VESPA** and **VESPAl**:

##### `.h5` files

Multiple script generate `.h5` files. These files follow the [`hdf5`-standard](https://www.hdfgroup.org/solutions/hdf5)
and can be processed in python using the library [`h5py`](https://www.h5py.org/). Generally the files are segmented into datasets that can be accessed using the protein accession in the fasta file. Each dataset is matrix shaped and usually has an `NxM`-shape, where `N` is the length of the respective protein length and `M` is the possible number of mutants (including self). `M` is a constant and the mutant length and mutant-order are determined by `MUTANT_ORDER` in `predict/config.py`. Empty fields that were not calculated/ specified contain a -1.

##### Mutations file

A simple text file with a protein ID and one mutation per-line (i.e., mutations separated by \n for newline).
Every mutation should be specified by `<PROTEIN_ID>_<SAV-String>` separated by an underscore.

The sequence ID needs to be equivalent to the one in the sequence FASTA file.
The SAV string has the format: `<Original Amino Acid><Position><Replacement Amino Acid>`

Example:

```txt
ENSP00000355206_I1L
ENSP00000355206_I1V
ENSP00000355206_I1L
ENSP00000355206_I1K
ENSP00000355206_I1T
```

##### Log odds ratio output

When running `--single_csv data/single_logodds.csv`, the log odds ratio output file contains the mutations determined by `<PROTEIN_ID>_<SAV-String>`, followed by a `;` and the log odds ratio for a particular mutation.

Example:

```csv
B3VI55_LIPSTSTABLE_A438Q;0.053829677402973175
B3VI55_LIPSTSTABLE_A438N;0.061238136142492294
B3VI55_LIPSTSTABLE_A438Y;0.012603843584656715
B3VI55_LIPSTSTABLE_A438M;0.012212143279612064
B3VI55_LIPSTSTABLE_A438H;0.0241163931787014
B3VI55_LIPSTSTABLE_A438W;0.004520446062088013
```

In case the precitions are written into a directory, e.g. by specifying `--csv_dir` in the `vespa_logodds`, the script will create one file per sequence, named by sequence ID (note: the ID will be normalized, i.e each special char will be replaced by `_`).
The file will contain one `<Mutation-String>;score` per-line.

Example file called `B3VI55_LIPSTSTABLE`:

```csv
A438Q;0.053829677402973175
A438N;0.061238136142492294
A438Y;0.012603843584656715
A438M;0.012212143279612064
A438H;0.0241163931787014
A438W;0.004520446062088013
```

##### VESPA and VESPAl output

**VESPA** and **VESPAl** will both generate one csv file per protein in the specified output directory (with `--output`).
To circumvent naming issues due to long sequence ID's, the csv files will be numbered by sequence occurence in the FASTA file. A lookup file `map.json` will be created in the output directory containing a dictionary mapping from number to sequence ID.

Example `map.json`:

```json
{
    "0": "B3VI55_LIPSTSTABLE",
    "1": "BF520_ENV",
    "2": "BG_STRSQ",
    "3": "BG505_ENV",
    "4": "HG_FLU",
    "5": "MTH3_HAEAESTABILIZED"
}
```

The individual files will contain rows with the mutations along with the respective predictions of `VESPA`, `VESPAl`, or both models.

Example `0.csv`:

```csv
Mutant;VESPAl;VESPA
M0A;0.4457732174287125;0.3520255108578212
M0L;0.3191178420567241;0.2717188481387661
M0G;0.5355136080284415;0.4110670843315182
M0V;0.3594337197937546;0.2971971641898669
M0S;0.4457732174287125;0.35202555423053
M0R;0.4457732174287125;0.35202621644931126
```

## Development Roadmap

- [ ] Write comprehensive tests
- [x] Publish pypi package
- [x] Install from github release
- [ ] Contributing

## Installation from current Github Release

**WARNING Experimental**: To install the current release from github you can use:

```bash
python -m pip install https://github.com/Rostlab/VESPA/releases/download/v0.9.0-beta/vespa-0.9.0b0.tar.gz
```

## Cite

If you want to credit us, feel free to cite

Marquet, C., Heinzinger, M., Olenyi, T. et al. Embeddings from protein language models predict conservation and variant effects. Hum Genet (2021). <https://doi.org/10.1007/s00439-021-02411-y>

```Bibtex
@article{Marquet2021,
  doi = {10.1007/s00439-021-02411-y},
  url = {https://doi.org/10.1007/s00439-021-02411-y},
  year = {2021},
  month = dec,
  publisher = {Springer Science and Business Media {LLC}},
  author = {C{\'{e}}line Marquet and Michael Heinzinger and Tobias Olenyi and Christian Dallago and Kyra Erckert and Michael Bernhofer and Dmitrii Nechaev and Burkhard Rost},
  title = {Embeddings from protein language models predict conservation and variant effects},
  journal = {Human Genetics}
}
```
