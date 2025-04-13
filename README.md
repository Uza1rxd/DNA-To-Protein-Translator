# DNA to Protein Translator

A user-friendly GUI application built with Python and Tkinter that translates nucleotide sequences (DNA/RNA) into protein sequences. The tool is designed for molecular biologists and bioinformaticians to easily explore reading frames, visualize codon breakdowns, and highlight motifs within translated proteins.

---

## Features

- **Sequence Input Options:**
  - Paste or type nucleotide sequences directly into the text box.
  - Load sequences from `.txt` or `.fasta` files (FASTA headers are automatically removed).

- **Input Sanitization & RNA Support:**
  - Automatically removes invalid characters and extraneous whitespace.
  - Detects RNA input (with U instead of T) and converts it for seamless translation.

- **Reading Frame Selection:**
  - Translate in a specific reading frame using radio buttons (+1, +2, +3, -1, -2, -3).
  - Option to translate all six reading frames simultaneously.

- **Codon-by-Codon Breakdown:**
  - Displays a detailed mapping of each codon (with its starting position) to its corresponding amino acid.

- **Genetic Code Options:**
  - Choose between the Standard genetic code or the Mitochondrial genetic code for translation.

- **Visual Enhancements:**
  - Highlights stop codons (`*`) in red.
  - Highlights motifs—regions starting with a start codon (M) and ending with a stop codon (`*`)—with a yellow background.
  - Highlights start codons (ATG/AUG) within the input sequence.

- **Result Options:**
  - Easily copy the translated sequence to the clipboard.
  - Export translated results to a `.txt` or `.fasta` file.

---

## Getting Started

### Prerequisites

- **Python 3.7+**
- **Tkinter** (usually included with Python)
- No additional third-party modules are required.

### Installation

1. **Clone the Repository:**

    ```bash
    git clone https://github.com/yourusername/dna-to-protein-translator.git
    cd dna-to-protein-translator
    ```

2. **(Optional) Create and Activate a Virtual Environment:**

    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows, use: venv\Scripts\activate
    ```

3. **Run the Application:**

    ```bash
    python translator_gui.py
    ```

---

## Usage Instructions

1. **Input Sequence:**
   - Paste or type your nucleotide sequence (DNA or RNA) in the provided text area.
   - Alternatively, load a sequence using the "Upload File (.txt/.fasta)" button.

2. **Select Options:**
   - Choose the desired reading frame using the radio buttons.
   - Check "Translate All 6 Frames" if you’d like to see translations for all frames.
   - Optionally, enable "Show reverse complement in input" to view the reverse complement.

3. **Genetic Code:**
   - Select "Standard" or "Mitochondrial" from the dropdown menu to set the appropriate genetic code.

4. **Translate:**
   - Click the "Translate Sequence" button.
   - The resulting protein sequence, along with a codon-by-codon breakdown, will appear in the output area.
   - The tool highlights key features such as stop codons in red and motifs (from ‘M’ to ‘*’) in yellow.

5. **Copy or Export Results:**
   - Use the "Copy Result" button to copy the translation to your clipboard.
   - Use the "Export Result" button to save the output as a `.txt` or `.fasta` file.

---

