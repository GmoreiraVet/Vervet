#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OCR GenBank Accession Extractor
"""

import re
import os
import logging
import cv2
import numpy as np
from easyocr import Reader
from Bio import Entrez, SeqIO
from itertools import product
from io import StringIO

# Set up logging. 
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
# Entrez requires an email address for identification. Replace with your actual email. (or don't worry about it)
Entrez.email = "your_email@example.com"

# Substitution dictionary for OCR error correction. Unfortunately, OCR can misread characters that look similar, especially in accession numbers. This dictionary maps commonly confused characters to their possible alternatives. In the future, I'd like to remove unnecessary substitutions (preventing dalse positives and reducing the number of variants to test) and add more based on observed OCR errors.
SUBSTITUTIONS = {
    'O': ['O', '0', 'Q'],
    '0': ['O', '0', 'Q'],
    'Q': ['O', '0', 'Q'],
    'I': ['I', '1'],
    '1': ['I', '1'],
    'S': ['S', '5'],
    '5': ['S', '5'],
    'B': ['B', '8'],
    '8': ['B', '8'],
    'A': ['A', '4'],
    '4': ['A', '4', 'K'],
    'T': ['T', '7'],
    '7': ['T', '7'],
    'G': ['G', '6'],
    '6': ['G', '6'],
    'Z': ['Z', '2'],
    '2': ['Z', '2'],
    'E': ['E', '3'],
    '3': ['E', '3'],
    'K': ['K', '4'],
}

# Looser pattern to match accession-like strings (not anchored). 2 letters followed by 5 or 6 digits, optionally with a version suffix. This will be used to generate candidate variants from OCR tokens.
ACCESSION_PATTERN = re.compile(r'[A-Z]{1,2}\d{5,6}(\.\d+)?', re.IGNORECASE)

def extract_text_with_boxes(image_path):
    reader = Reader(['en'], gpu=False) # Here Gpu can be turned on, to be changed in future versions to be more dynamic (detect if gpu is available or not).
    logging.info("Extracting text from the image...")
    results = reader.readtext(image_path, detail=1) #detail=1 returns bounding box, text, and confidence score for each detected text element. This allows us to draw boxes around the detected text in the image later on.
    return results

# This function takes the original image and the OCR results (which include bounding boxes and detected text) and draws green boxes around the detected text. It also annotates the image with the detected text for easier manual review. The resulting image is saved as "ocr_highlighted.png". For debugging, to be removed.
def draw_ocr_boxes(image_path, ocr_results, output_path="ocr_highlighted.png"):
    image = cv2.imread(image_path)
    for (bbox, text, _) in ocr_results:
        pts = np.array(bbox).astype(int).reshape((-1, 1, 2))
        cv2.polylines(image, [pts], isClosed=True, color=(0, 255, 0), thickness=2)
        cv2.putText(image, text, (pts[0][0][0], pts[0][0][1] - 10),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 255, 0), 2)
    cv2.imwrite(output_path, image)
    logging.info(f"OCR-highlighted image saved as: {output_path}")

# This function takes the OCR results, extracts just the detected text, and saves it to a plain text file (default "ocr_text_output.txt"). It also returns a list of the extracted text lines for further processing. This allows us to have a clean text file with just the OCR results, which can be useful for debugging and manual review. 
def save_raw_text(ocr_results, output_txt="ocr_text_output.txt"):
    lines = [text for _, text, _ in ocr_results]
    with open(output_txt, 'w') as f:
        f.write('\n'.join(lines))
    logging.info(f"OCR text saved to: {output_txt}")
    return lines

# This function generates all possible variants of a given word by substituting characters based on the SUBSTITUTIONS dictionary. It first converts the input word to uppercase, then creates a list of possible substitutions for each character. The itertools.product function is used to generate all combinations of these substitutions. Finally, it filters the generated variants to include only those that match the ACCESSION_PATTERN, ensuring that we only test valid accession-like strings against the Entrez database.
def generate_variants(word):
    word = word.upper()
    options = [SUBSTITUTIONS.get(c, [c]) for c in word]
    variants = [''.join(p) for p in product(*options)]
    return [v for v in variants if ACCESSION_PATTERN.fullmatch(v)]

# This function takes an accession number, queries the NCBI Entrez database to check if it is valid, and returns the FASTA content if it is. If the accession number is not valid or if there is an error during the fetch, it returns None. This function is crucial for validating the candidate accession numbers generated from the OCR text.
def is_accession_valid(acc):
    try:
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
        content = handle.read()
        handle.close()
        return content if content.startswith('>') else None
    except Exception as e:
        logging.debug(f"Entrez fetch failed for {acc}: {e}")
        return None

# This function takes a line of text, removes most non-alphanumeric characters (except for dots, which are used in accession numbers for versioning), and splits the cleaned line into tokens. This helps us extract potential accession number candidates from the OCR text while ignoring irrelevant characters.
def tokenize_line(line):
    # Strip most non-alphanum (keep dots for versioning)
    cleaned = re.sub(r'[^A-Za-z0-9.]', ' ', line)
    return cleaned.split()

# This function takes a list of text lines, combines them into a single string, and uses a regular expression to extract unique accession number candidates. The regex looks for patterns that match typical GenBank accession numbers (1-2 letters followed by 5-6 digits, optionally with a version suffix). This provides an initial set of candidates that we can then validate against the Entrez database.
def extract_accession_candidates(text_lines):
    full_text = " ".join(text_lines)
    return list(set(re.findall(r'[A-Z]{1,2}\d{5,6}(?:\.\d+)?', full_text, re.IGNORECASE)))

# checks each candidate against the Entrez database using the is_accession_valid function. Valid accessions are stored in a dictionary, while unresolved tokens (those that could not be validated) are stored separately for potential manual review. The function also keeps track of the total number of tokens processed for reporting purposes.
def find_valid_accessions(text_lines):
    logging.info("Generating and validating accession numbers...")
    valid_accessions = {}
    unresolved_tokens = {}
    tested = set()
    total_tokens = 0

    for line in text_lines:
        tokens = tokenize_line(line)
        total_tokens += len(tokens)
        for token in tokens:
            upper_token = token.upper()
            if upper_token in valid_accessions or upper_token in unresolved_tokens:
                continue

            token_no_ver = upper_token.split('.')[0]
            if len(token_no_ver) not in (7, 8):
                continue

            variants = generate_variants(token_no_ver)

            if ACCESSION_PATTERN.fullmatch(token_no_ver):
                candidates = [token_no_ver] + [v for v in variants if v != token_no_ver]
            else:
                candidates = variants

            if not candidates:
                continue

            found = False
            for candidate in candidates:
                if candidate in tested:
                    continue
                tested.add(candidate)
                fasta = is_accession_valid(candidate)
                if fasta:
                    valid_accessions[upper_token] = (candidate, fasta)
                    logging.info(f"Valid accession: {upper_token} → {candidate}")
                    found = True
                    break
            if not found:
                unresolved_tokens[upper_token] = variants
    return valid_accessions, unresolved_tokens, total_tokens

# This function allows the user to manually review unresolved tokens that could not be automatically validated. For each unresolved token, it displays the token and up to 5 suggested variants based on the OCR error correction. The user can then input a corrected accession number, which is validated against the Entrez database. If valid, it is added to the list of valid accessions; if not, it is skipped. This manual review process provides an opportunity to recover valid accession numbers that may have been missed due to OCR errors or limitations in the automatic validation process.
def manual_review(unresolved_tokens, valid_accessions):
    logging.info("Manual review of unresolved accession numbers:")
    for token, variants in unresolved_tokens.items():
        if not variants:
            continue
        print(f"\nUnresolved token: {token}")
        print("Suggested variants (up to 5):")
        print(", ".join(variants[:5]))

        manual_input = input("Enter correct accession (or press Enter to skip): ").strip().upper()
        if manual_input:
            if not ACCESSION_PATTERN.fullmatch(manual_input):
                print(f"⚠️ '{manual_input}' does not match accession format, skipping.")
                continue
            fasta = is_accession_valid(manual_input)
            if fasta:
                valid_accessions[token] = (manual_input, fasta)
                logging.info(f"Manually added valid accession: {token} → {manual_input}")
            else:
                print(f"⚠️ Accession '{manual_input}' not found in Entrez. Skipping.")
                logging.warning(f"Manual accession fetch failed: {manual_input}")
        else:
            logging.info(f"Skipped manual input for token: {token}")

# This function takes the dictionary of valid accessions (where the key is the original OCR token and the value is a tuple of the validated accession and its FASTA content) and saves all the valid sequences to a specified FASTA file. It uses Biopython's SeqIO to write the sequences in FASTA format. If there are no valid accessions, it logs a warning and does not create an output file.
def fetch_and_save_sequences(accession_dict, output_fasta_path):
    if not accession_dict:
        logging.warning("No valid accession numbers to fetch.")
        return
    records = []
    for orig, (acc, fasta_text) in accession_dict.items():
        handle = StringIO(fasta_text)
        seq = next(SeqIO.parse(handle, "fasta"), None)
        if seq:
            records.append(seq)
    SeqIO.write(records, output_fasta_path, "fasta")
    logging.info(f"Saved {len(records)} sequences to {output_fasta_path}")


# The main function orchestrates the entire process. It prompts the user for the image file path and the output FASTA file path. It then checks if the image file exists, extracts text using OCR, draws bounding boxes on the image for visualization, saves the raw OCR text to a file, and processes the text to find valid accession numbers. It reports the number of tokens detected and valid accessions found, and offers the user an option to manually review unresolved tokens. Finally, it fetches the valid sequences and saves them to a FASTA file.
def main():
    image_path = input("Enter path to image file (e.g., image.jpg): ").strip()
    output_fasta_path = input("Enter output FASTA file path (e.g., output.fasta): ").strip()

    if not os.path.exists(image_path):
        logging.error(f"Image file not found: {image_path}")
        return

    ocr_results = extract_text_with_boxes(image_path)
    draw_ocr_boxes(image_path, ocr_results)
    lines = save_raw_text(ocr_results)

    valid_accessions, unresolved_tokens, total_tokens = find_valid_accessions(lines)

    print(f"\nDetected {total_tokens} token(s) from OCR text.")
    print(f"Found {len(valid_accessions)} valid accession number(s) out of these tokens.")

    logging.info(f"✅ Vervet successfully identified {len(valid_accessions)} valid accession number(s).")
    logging.info(f"🔎 {len(unresolved_tokens)} token(s) could not be automatically resolved.")

    if unresolved_tokens:
        response = input(
            f"\n{len(unresolved_tokens)} potential accession number candidates could not be automatically resolved.\n"
            "Would you like to review and attempt to correct them manually? (yes/no): "
        ).strip().lower()
        if response in ("yes", "y"):
            manual_review(unresolved_tokens, valid_accessions)
        else:
            logging.info("Manual review skipped by user.")

    fetch_and_save_sequences(valid_accessions, output_fasta_path)

    print("Done!")

if __name__ == "__main__":
    main()
