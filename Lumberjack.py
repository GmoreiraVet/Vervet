import easyocr
import re
from Bio import Entrez

def extract_text_from_image(image_path):
    """Extract text from the image using EasyOCR."""
    try:
        reader = easyocr.Reader(['en'])  # Specify English language
        result = reader.readtext(image_path)
        
        # Extract text content from OCR result
        extracted_text = " ".join([res[1] for res in result])
        return extracted_text
    except Exception as e:
        print(f"Error extracting text from image: {str(e)}")
        return None

def adjust_and_validate_accession_numbers(text):
    """Correct misread characters and validate potential accession numbers."""
    # Extract possible accession numbers (1 letter + 5 digits OR 2 letters + 6 digits)
    potential_accessions = re.findall(r'\b[a-zA-Z]{1,2}\d{5,6}\b', text)
    adjusted_accessions = []
    valid_accessions = []
    
    for acc in potential_accessions:
        # Correct possible misinterpretations
        first_part = re.sub(r'0', 'O', acc[:2], flags=re.IGNORECASE)  # Convert '0' to 'O' in first two characters
        second_part = re.sub(r'O', '0', acc[2:], flags=re.IGNORECASE)  # Convert 'O' to '0' in the rest
        
        adjusted_acc = first_part + second_part
        adjusted_accessions.append(adjusted_acc)
        
        # Validate format: either 1 letter + 5 digits OR 2 letters + 6 digits
        if re.match(r'^[A-Za-z]{1}\d{5}$', adjusted_acc) or re.match(r'^[A-Za-z]{2}\d{6}$', adjusted_acc):
            valid_accessions.append(adjusted_acc)
    
    return adjusted_accessions, valid_accessions

def validate_accession_number(accession_number):
    """Check if the accession number exists in GenBank."""
    try:
        handle = Entrez.esearch(db="nucleotide", term=accession_number, retmax=1)
        record = Entrez.read(handle)
        handle.close()
        return record['Count'] != '0'
    except Exception as e:
        print(f"Error validating accession number {accession_number}: {str(e)}")
        return False

def fetch_fasta_from_genbank(accession_numbers):
    """Fetch FASTA sequences from GenBank for the provided accession numbers."""
    Entrez.email = "your.email@example.com"  # Replace with your email
    
    sequences = []
    failed_accessions = []
    
    for acc_num in accession_numbers:
        if not validate_accession_number(acc_num):
            print(f"Invalid accession number: {acc_num}. Skipping.")
            failed_accessions.append(acc_num)
            continue
        
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc_num, rettype="fasta", retmode="text")
            fasta_data = handle.read()
            sequences.append(fasta_data)
            handle.close()
        except Exception as e:
            print(f"Failed to fetch {acc_num}: {str(e)}")
            failed_accessions.append(acc_num)
    
    return sequences, failed_accessions

def save_fasta(sequences, output_file):
    """Save the fetched sequences to a FASTA file."""
    try:
        with open(output_file, 'w') as f:
            for seq in sequences:
                f.write(seq + "\n")
        print(f"FASTA file saved to {output_file}.")
    except Exception as e:
        print(f"Error saving FASTA file: {str(e)}")

def main():
    image_path = input("Enter the path to the image file (e.g., image.jpg, image.png): ")
    output_fasta_path = input("Enter the output FASTA file path (e.g., output.fasta): ")

    print("Extracting text from the image...")
    text = extract_text_from_image(image_path)
    
    if not text:
        print("Failed to extract text from the image.")
        return

    print("Adjusting and validating accession numbers...")
    _, valid_accessions = adjust_and_validate_accession_numbers(text)
    
    if not valid_accessions:
        print("No valid accession numbers found in the image.")
        return

    print(f"Detected {len(valid_accessions)} accession numbers: {valid_accessions}")
    print("Fetching sequences from GenBank...")
    sequences, failed_accessions = fetch_fasta_from_genbank(valid_accessions)
    
    if sequences:
        save_fasta(sequences, output_fasta_path)
    else:
        print("No sequences were successfully fetched.")
    
    if failed_accessions:
        print(f"Failed to fetch the following accession numbers: {failed_accessions}")

if __name__ == "__main__":
    main()
