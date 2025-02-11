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

def extract_accession_numbers(text):
    """Extract valid accession numbers from text using regex."""
    # Regex pattern to match GenBank accession number formats (RefSeq, Genomic, mRNA, etc.)
    accession_pattern = r'\b([A-Za-z]{2,4}_?\d{6,9}[A-Za-z]?)\b'  # This matches most common GenBank formats
    return re.findall(accession_pattern, text)

def validate_accession_number(accession_number):
    """Check if the accession number exists in GenBank."""
    try:
        handle = Entrez.esearch(db="nucleotide", term=accession_number, retmax=1)
        record = Entrez.read(handle)
        handle.close()
        # If no results are returned, it's an invalid accession number
        if record['Count'] == '0':
            return False
        return True
    except Exception as e:
        print(f"Error validating accession number {accession_number}: {str(e)}")
        return False

def fetch_fasta_from_genbank(accession_numbers):
    """Fetch FASTA sequences from GenBank for the provided accession numbers."""
    Entrez.email = "your.email@example.com"  # Replace with your email
    
    sequences = []
    failed_accessions = []
    
    for acc_num in accession_numbers:
        # Validate the accession number before fetching data
        if not validate_accession_number(acc_num):
            print(f"Invalid accession number: {acc_num}. Skipping.")
            failed_accessions.append(acc_num)
            continue
        
        try:
            # Query GenBank for the accession number
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
    # Input paths
    image_path = input("Enter the path to the image file (e.g., image.jpg, image.png): ")
    output_fasta_path = input("Enter the output FASTA file path (e.g., output.fasta): ")

    # Extract text from the image
    print("Extracting text from the image...")
    text = extract_text_from_image(image_path)
    
    if not text:
        print("Failed to extract text from the image.")
        return

    # Extract accession numbers from the text
    print("Extracting accession numbers from the text...")
    accession_numbers = extract_accession_numbers(text)

    if not accession_numbers:
        print("No accession numbers found in the image.")
        return

    print(f"Found accession numbers: {accession_numbers}")

    # Fetch sequences from GenBank
    print("Fetching sequences from GenBank...")
    sequences, failed_accessions = fetch_fasta_from_genbank(accession_numbers)

    if sequences:
        # Save the fetched sequences to a FASTA file
        save_fasta(sequences, output_fasta_path)
    else:
        print("No sequences were successfully fetched.")

    # Report failed accessions
    if failed_accessions:
        print(f"Failed to fetch the following accession numbers: {failed_accessions}")

if __name__ == "__main__":
    main()
