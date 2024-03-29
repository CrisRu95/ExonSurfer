
# ExonSurfer

<p align="center">
  <img src="./ExonSurfer/resources/ES.png" width="100" height="100">
</p>

ExonSurfer is a python3 program, designed for the design of high specific transcript primers, using the Ensembl database and the blastn algorithm in order to identify and find the proper exon binding of the transcripts.

The advantages of ExonPrimerSurfer are:

* Fast and accurate primer design, with user specified parameters like primer size, GC content, etc
* User friendly [WEB](https://exonsurfer.i-med.ac.at/) to setup search options and visualize results
* Highly specific primers, with low non-specific binding
* Fully automated workflows

The program is easy to install, simply clone the repository to your local machine and install the required dependencies. 
## Installation

Clone the repository:

```git
git clone https://github.com/CrisRu95/ExonSurfer.git

```

Install the dependencies:

```bash
pip install -r requirements.txt
```

Install the package

```bash
python3 setup.py install --user
```

## Usage

Once the program is installed, it can be used to design primers for any given gene transcript.

`exon_surfer.py --gene <gene_name> --transcript <transcript_name> --out <path_out> --release <release> --design_dict <design_dict> --files <save_files> --e_value <e_value> --i_cutoff <i_cutoff> --max_sep <max_sep>
`

This will output a list of primers that can be used for the specified gene transcript.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
