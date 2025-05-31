# Image to GenBank FASTA Downloader
![plot](Untitled-1.png)

Vervet extracts accession numbers from an image using Optical Character Recognition (OCR) with EasyOCR, validates them, and then fetches the corresponding sequences from GenBank. It saves the sequences to a FASTA file.

## Features:
- Extracts text from an image (e.g., `.jpg`, `.png`) containing accession numbers.
- Uses EasyOCR to perform OCR and extract accession numbers.
- Validates accession numbers to check if they exist in GenBank.
- Substitutes commonly misinterpreted characters when a candidate accession number doesn´t exist in genbank
- Fetches the corresponding sequences from GenBank using the Biopython library.
- Saves the fetched sequences in a FASTA format file.

## Dependencies

1. **Python 3.x**: Ensure you have Python 3 installed on your machine.
2. **EasyOCR**: For extracting text (accession numbers) from images using OCR.
3. **Biopython**: For interacting with GenBank and fetching sequences.
4. **Requests**: For handling HTTP requests to GenBank.
5. **Pillow**: For image processing.
   
You can install these dependencies by running the following command:

```bash
pip install easyocr biopython requests pillow
```

# Citation
If you use this code in your work, please cite it using the following information:
```
cff-version: 1.2.0
message: "If you use this code in your work, please cite it using the following information."
authors:
  - name: "Guilherme Moreira"
    orcid: "https://orcid.org/0009-0002-2828-7202"
    affiliation: "ICBAS - Instituto de Ciências Biomédicas Abel Salazar"
title: "Lumberjack"
version: "1.0.0"
doi: "https://doi.org/10.5281/zenodo.14850245"
date-released: 2025-02-11
repository: "https://github.com/GmoreiraVet/Lumberjack"
```
The DOI for citing this work is: 10.5281/zenodo.14850245.
