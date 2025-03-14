Sure! Here’s the updated README file with the requested modification:

# MTase Classification Pipeline

This repository provides a pipeline for the classification of Methyltransferases (MTases) using Hidden Markov Models (HMMs). The pipeline consists of three main steps: installing the required packages, conducting HMMer searches, and detecting regions for classification.

## Table of Contents

1. [Installation](#installation)
2. [HMMer Search](#hmmer-search)
3. [Region Detection](#region-detection)
4. [Classification](#classification)
5. [Web Application](#web-application)
6. [Additional Information](#additional-information)

## Installation

Before running the pipeline, you will need to install the required packages. Execute the following commands in your Colab Notebook:

```bash
!git clone https://github.com/MVolobueva/MTase-classification.git
!sudo apt-get -y install hmmer
!git clone https://github.com/isrusin/etsv
!python3 -m pip install -e etsv
```

The visualization of the steps in the pipeline is shown in the image below:

```python
from IPython.display import Image
Image('/content/MTase-classification/Pipeline.jpg')
```

## HMMer Search

Conduct the HMMer search to identify Methyltransferase sequences. Use the following command:

```bash
!hmmsearch --cpu 3 -E 0.01 --domE 0.01 --incE 0.01 --incdomE 0.01 \
        -o /dev/null --noali -A file.stk \
        /content/MTase-classification/HMM_profiles/selected_profiles.hmm /content/MTase-classification/Sample_MTases/MTase_sequences.fasta
```

## Region Detection

After running the HMMer search, the next step is to detect regions in the alignment. Use this command:

```bash
!./MTase-classification/Scripts/get_aln_regions.py \
  /content/MTase-classification/profile-markup/All_profile_region.csv \
  /content/file.stk > region_alignments.tsv
```

## Classification

Finally, perform the classification of the detected regions with the following command:

```bash
!python ./MTase-classification/Scripts/classification.py \
  --t /content/region_alignments.tsv \
  --m several_cat_domains.tsv \
  --c class.tsv
```

## Web Application

To illustrate the workings of the pipeline, we have developed a web application. You can access it at the following link:

[MTase Classification Web Application](https://mtase-pipeline-6g1yfq9ugw8.streamlit.app/)

## Additional Information

For detailed information about the classification method used in this pipeline, please refer to the latest version of the manuscript available at: [Classification Manuscript](https://www.biorxiv.org/content/10.1101/2023.12.13.571470v4.full).

## Contributing

Contributions to this project are welcome! Please feel free to fork the repository and submit pull requests.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

We would like to acknowledge the developers and contributors of the libraries and tools used in this project. Special thanks to the HMMer team for their contribution to bioinformatics.

---

Feel free to modify any sections as per your needs!
