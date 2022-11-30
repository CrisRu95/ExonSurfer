
# ExonPrimerSurfer


ExonPrimerSurfer is a python3 program, designed for the design of high specific transcript primers, using the Ensembl database and the blastn algorithm in order to identify and find the proper exon binding of the transcripts.

The advantages of ExonPrimerSurfer are:

* Fast and accurate primer design, with user specified parameters like primer size, GC content, etc
* User friendly GUI to setup search options and visualize results
* Highly specific primers, with low non-specific binding
* Fully automated workflows

The program is easy to install, simply clone the repository to your local machine and install the required dependencies. 
## Installation

Clone the repository:

```git
git clone https://github.com/username/ExonPrimerSurfer.git
```

Install the dependencies:

```bash
pip install -r requirements.txt
```

## Usage

Once the program is installed, it can be used to design primers for any given gene transcript.

`python exonprimersurfer.py --gene <gene_name> --transcript <transcript_name>`

This will output a list of primers that can be used for the specified gene transcript.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
