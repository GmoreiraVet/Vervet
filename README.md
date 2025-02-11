# Image to GenBank FASTA Downloader

This script extracts accession numbers from an image using Optical Character Recognition (OCR) with EasyOCR, validates them, and then fetches the corresponding sequences from GenBank. It saves the sequences to a FASTA file.

## Features:
- Extracts text from an image (e.g., `.jpg`, `.png`) containing accession numbers.
- Uses EasyOCR to perform OCR and extract accession numbers.
- Validates accession numbers to check if they exist in GenBank.
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
