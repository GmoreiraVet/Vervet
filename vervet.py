import re
import os
import logging
import cv2
import numpy as np
from easyocr import Reader
from Bio import Entrez, SeqIO
from itertools import product
from io import StringIO

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
Entrez.email = "your_email@example.com"  # Replace with your actual email

# Substitution dictionary for OCR error correction
SUBSTITUTIONS = {
    'O': ['O', '0'],
    '0': ['O', '0'],
    'I': ['I', '1'],
    '1': ['I', '1'],
    'S': ['S', '5'],
    '5': ['S', '5'],
    'B': ['B', '8'],
    '8': ['B', '8'],
    'A': ['A', '4'],
    '4': ['A', '4'],
    'T': ['T', '7'],
    '7': ['T', '7'],
    'G': ['G', '6'],
    '6': ['G', '6'],
    'Z': ['Z', '2'],
    '2': ['Z', '2'],
    'E': ['E', '3'],
    '3': ['E', '3'],
    'K': ['K', '4'],
    '4': ['4', 'K'],
}

ACCESSION_PATTERN = re.compile(r'^[A-Z]{1,2}\d{5,6}(\.\d+)?$')
CANDIDATE_PATTERN = re.compile(r'^[A-Z]{1,2}[A-Z0-9]{4,7}\d$')

VALID_START_CHARS = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789O0I1S5B8A4T7G6Z2E3K4')

def extract_text_with_boxes(image_path):
    reader = Reader(['en'], gpu=False)
    logging.info("Extracting text from the image...")
    results = reader.readtext(image_path, detail=1)
    return results

def draw_ocr_boxes(image_path, ocr_results, output_path="ocr_highlighted.png"):
    image = cv2.imread(image_path)
    for (bbox, text, _) in ocr_results:
        pts = np.array(bbox).astype(int).reshape((-1, 1, 2))
        cv2.polylines(image, [pts], isClosed=True, color=(0, 255, 0), thickness=2)
        cv2.putText(image, text, (pts[0][0][0], pts[0][0][1] - 10),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 255, 0), 2)
    cv2.imwrite(output_path, image)
    logging.info(f"OCR-highlighted image saved as: {output_path}")

def save_raw_text(ocr_results, output_txt="ocr_text_output.txt"):
    lines = [text for _, text, _ in ocr_results]
    with open(output_txt, 'w') as f:
        f.write('\n'.join(lines))
    logging.info(f"OCR text saved to: {output_txt}")
    return lines

def generate_variants(word):
    word = word.upper()
    options = [SUBSTITUTIONS.get(c, [c]) for c in word]
    variants = [''.join(p) for p in product(*options)]
    return [v for v in variants if ACCESSION_PATTERN.fullmatch(v)]

def is_accession_valid(acc):
    try:
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
        content = handle.read()
        handle.close()
        return content if content.startswith('>') else None
    except Exception as e:
        logging.debug(f"Entrez fetch failed for {acc}: {e}")
        return None

def find_valid_accessions(text_lines):
    logging.info("Generating and validating accession numbers...")
    valid_accessions = {}
    unresolved_tokens = {}
    tested = set()

    for line in text_lines:
        tokens = re.findall(r'\w+', line)
        for token in tokens:
            upper_token = token.upper()
            if upper_token in valid_accessions or upper_token in unresolved_tokens:
                continue
            if len(upper_token) < 6 or not re.search(r'[A-Z]', upper_token):
                continue
            variants = generate_variants(upper_token)
            found = False
            for variant in variants:
                if variant in tested:
                    continue
                tested.add(variant)
                fasta = is_accession_valid(variant)
                if fasta:
                    valid_accessions[upper_token] = (variant, fasta)
                    logging.info(f"Valid accession: {upper_token} → {variant}")
                    found = True
                    break
            if not found:
                unresolved_tokens[upper_token] = variants

    # Filter unresolved tokens: only keep those that
    # - match candidate pattern,
    # - start with allowed chars,
    # - and have at least one variant suggestion
    unresolved_tokens = {
        tok: vars for tok, vars in unresolved_tokens.items()
        if CANDIDATE_PATTERN.fullmatch(tok)
        and tok[0] in VALID_START_CHARS
        and len(vars) > 0
    }

    return valid_accessions, unresolved_tokens

def manual_review(unresolved_tokens, valid_accessions):
    logging.info("Manual review of unresolved accession numbers:")
    for token, variants in unresolved_tokens.items():
        if not CANDIDATE_PATTERN.fullmatch(token):
            continue
        print(f"\nUnresolved token: {token}")
        print("Suggested variants (up to 5):")
        print(", ".join(variants[:5]) if variants else "None")

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

def main():
    image_path = input("Enter path to image file (e.g., image.jpg): ").strip()
    output_fasta_path = input("Enter output FASTA file path (e.g., output.fasta): ").strip()

    if not os.path.exists(image_path):
        logging.error(f"Image file not found: {image_path}")
        return

    ocr_results = extract_text_with_boxes(image_path)
    draw_ocr_boxes(image_path, ocr_results)
    lines = save_raw_text(ocr_results)

    valid_accessions, unresolved_tokens = find_valid_accessions(lines)

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

if __name__ == "__main__":
    main()
