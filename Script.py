import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
import os
import re

# Genetic code dictionaries (extended as needed)
STANDARD_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

MITOCHONDRIAL_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'M',  # ATA codes for M in vertebrate mitochondria
    'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': 'W',  # TGA codes for W in mitochondria
    'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': '*', 'AGG': '*',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def sanitize_sequence(seq):
    """
    Remove FASTA headers and invalid whitespace.
    Keep only valid nucleotide characters: A, T, C, G, U.
    Additionally, if the sequence is RNA (contains U but not T), convert U to T.
    """
    lines = seq.splitlines()
    seq_clean = ""
    for line in lines:
        if line.startswith('>'):  # Skip FASTA header lines
            continue
        seq_clean += line.strip()
    seq_clean = seq_clean.upper().replace(" ", "").replace("\n", "")
    valid = set("ATCGU")
    seq_filtered = ''.join([ch for ch in seq_clean if ch in valid])
    # If the sequence seems to be RNA (contains U but no T), convert U to T for translation
    if "U" in seq_filtered and "T" not in seq_filtered:
        seq_filtered = seq_filtered.replace("U", "T")
    return seq_filtered


def reverse_complement(seq):
    """
    Compute the reverse complement of a nucleotide sequence.
    Handles both DNA (T) and RNA (U). If 'U' is in the sequence without 'T',
    it is treated as RNA.
    """
    if 'U' in seq and 'T' not in seq:
        complement = str.maketrans("AUGC", "UACG")
    else:
        complement = str.maketrans("ATCG", "TAGC")
    return seq.translate(complement)[::-1]


def translate_and_breakdown(seq, frame, code=STANDARD_CODE):
    """
    Translate a nucleotide sequence into a protein sequence for a given reading frame.
    Also generate a codon-by-codon breakdown.

    Returns a tuple (protein, breakdown_str).
    """
    protein = ""
    breakdown_lines = []
    if frame < 0:
        seq_mod = reverse_complement(seq)
        frame = abs(frame)
    else:
        seq_mod = seq
    seq_mod = seq_mod[frame - 1:]
    # Process each codon triplet
    for i in range(0, len(seq_mod) - 2, 3):
        codon = seq_mod[i:i + 3]
        if len(codon) < 3:
            break
        amino = code.get(codon, 'X')
        protein += amino
        # Display the codon number (starting position relative to original mod seq)
        breakdown_lines.append(f"{i + frame:3d}: {codon} -> {amino}")
    breakdown_str = "\n".join(breakdown_lines)
    return protein, breakdown_str


class TranslatorGUI:
    def __init__(self, master):
        self.master = master
        master.title("Nucleotide to Protein Translator")

        # Create frames for layout
        self.input_frame = tk.Frame(master)
        self.input_frame.pack(pady=5)

        self.options_frame = tk.Frame(master)
        self.options_frame.pack(pady=5)

        self.result_frame = tk.Frame(master)
        self.result_frame.pack(pady=5, fill="both", expand=True)

        # Input Section
        tk.Label(self.input_frame, text="Enter DNA/RNA Sequence:").grid(row=0, column=0, sticky="w")
        self.seq_text = scrolledtext.ScrolledText(self.input_frame, height=10, width=60)
        self.seq_text.grid(row=1, column=0, columnspan=2, padx=5, pady=5)

        self.upload_button = tk.Button(self.input_frame, text="Upload File (.txt/.fasta)", command=self.load_file)
        self.upload_button.grid(row=0, column=1, padx=5)

        # Options: Reading Frame Selection
        tk.Label(self.options_frame, text="Select Reading Frame:").grid(row=0, column=0, sticky="w")
        self.frame_var = tk.IntVar(value=1)
        frames = [("+1", 1), ("+2", 2), ("+3", 3), ("-1", -1), ("-2", -2), ("-3", -3)]
        col = 0
        for text, val in frames:
            rb = tk.Radiobutton(self.options_frame, text=text, variable=self.frame_var, value=val)
            rb.grid(row=1, column=col, sticky="w")
            col += 1

        # Option: Translate All 6 Frames
        self.all_frames_var = tk.BooleanVar(value=False)
        self.all_frames_check = tk.Checkbutton(self.options_frame, text="Translate All 6 Frames",
                                               variable=self.all_frames_var)
        self.all_frames_check.grid(row=2, column=0, columnspan=2, sticky="w")

        # Option: Automatically display reverse complement in input if needed
        self.rev_comp_var = tk.BooleanVar(value=False)
        self.rev_comp_check = tk.Checkbutton(self.options_frame, text="Show reverse complement in input",
                                             variable=self.rev_comp_var, command=self.highlight_input)
        self.rev_comp_check.grid(row=3, column=0, columnspan=2, sticky="w")

        # Genetic Code Selection
        tk.Label(self.options_frame, text="Genetic Code:").grid(row=4, column=0, sticky="w")
        self.code_var = tk.StringVar(value="Standard")
        self.code_option = tk.OptionMenu(self.options_frame, self.code_var, "Standard", "Mitochondrial")
        self.code_option.grid(row=4, column=1, sticky="w")

        # Translation button
        self.translate_button = tk.Button(self.options_frame, text="Translate Sequence", command=self.translate)
        self.translate_button.grid(row=4, column=2, padx=5, pady=5)

        # Output Section
        tk.Label(self.result_frame, text="Translation Results:").pack(anchor="w")
        self.result_text = scrolledtext.ScrolledText(self.result_frame, height=20, width=100)
        self.result_text.pack(padx=5, pady=5, fill="both", expand=True)

        # Copy and Export Buttons
        self.copy_button = tk.Button(master, text="Copy Result", command=self.copy_result)
        self.copy_button.pack(side="left", padx=10, pady=5)
        self.export_button = tk.Button(master, text="Export Result", command=self.export_result)
        self.export_button.pack(side="right", padx=10, pady=5)

        # Tag configuration for highlighting stop codons and motifs in output
        self.result_text.tag_config("stop", foreground="red")
        self.result_text.tag_config("motif", background="yellow")
        self.seq_text.tag_config("start", foreground="green")

    def load_file(self):
        file_path = filedialog.askopenfilename(title="Select Sequence File",
                                               filetypes=(("Text files", "*.txt *.fasta"), ("All files", "*.*")))
        if file_path:
            try:
                with open(file_path, "r") as file:
                    content = file.read()
                self.seq_text.delete("1.0", tk.END)
                self.seq_text.insert(tk.END, content)
                self.highlight_input()
            except Exception as e:
                messagebox.showerror("Error", f"Could not read file: {e}")

    def highlight_input(self):
        """
        Highlight start codons (ATG/AUG) in the input sequence.
        """
        self.seq_text.tag_remove("start", "1.0", tk.END)
        seq = self.seq_text.get("1.0", tk.END).upper()
        for codon in ["ATG", "AUG"]:
            start_index = "1.0"
            while True:
                pos = self.seq_text.search(codon, start_index, stopindex=tk.END)
                if not pos:
                    break
                end = f"{pos}+{len(codon)}c"
                self.seq_text.tag_add("start", pos, end)
                start_index = end

    def translate(self):
        """
        Sanitize the input, validate the sequence, translate using the chosen options,
        and display the resulting protein sequence(s) with codon-by-codon breakdowns.
        Also highlights motif regions (from a start 'M' to the next stop '*').
        """
        raw_seq = self.seq_text.get("1.0", tk.END)
        seq = sanitize_sequence(raw_seq)
        if len(seq) < 3:
            messagebox.showerror("Invalid Sequence", "The sequence is too short or contains no valid nucleotides.")
            return

        code = STANDARD_CODE if self.code_var.get() == "Standard" else MITOCHONDRIAL_CODE
        self.result_text.delete("1.0", tk.END)

        # Choose to translate a single chosen frame or all six frames
        if self.all_frames_var.get():
            frames = [1, 2, 3, -1, -2, -3]
        else:
            frames = [self.frame_var.get()]

        # Iterate over frames and generate output
        for frm in frames:
            protein, breakdown = translate_and_breakdown(seq, frm, code)
            frame_header = f"Frame {frm:+d}:\n"
            self.result_text.insert(tk.END, frame_header)
            self.result_text.insert(tk.END, protein + "\n")
            self.result_text.insert(tk.END, "Codon-by-Codon Breakdown:\n")
            self.result_text.insert(tk.END, breakdown + "\n")
            self.result_text.insert(tk.END, "-" * 60 + "\n\n")

        self.highlight_output()

    def highlight_output(self):
        """
        Highlight stop codons (*) in the protein output and motif regions
        (patches starting with 'M' and ending with '*' in the protein sequence).
        """
        # First remove previous tags
        self.result_text.tag_remove("stop", "1.0", tk.END)
        self.result_text.tag_remove("motif", "1.0", tk.END)

        # Highlight each stop codon individually
        index = "1.0"
        while True:
            pos = self.result_text.search(r"\*", index, stopindex=tk.END, regexp=True)
            if not pos:
                break
            end = f"{pos}+1c"
            self.result_text.tag_add("stop", pos, end)
            index = end

        # For motif highlighting, use regex on the entire text.
        # A motif is defined as a substring starting with 'M' (start) and ending with '*' (stop)
        full_text = self.result_text.get("1.0", tk.END)
        pattern = re.compile(r"M[^*]*\*")
        for match in pattern.finditer(full_text):
            start_idx = f"1.0+{match.start()}c"
            end_idx = f"1.0+{match.end()}c"
            self.result_text.tag_add("motif", start_idx, end_idx)

    def copy_result(self):
        result = self.result_text.get("1.0", tk.END)
        self.master.clipboard_clear()
        self.master.clipboard_append(result)
        messagebox.showinfo("Copied", "Result copied to clipboard!")

    def export_result(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".txt",
                                                 filetypes=[("Text files", "*.txt"), ("FASTA files", "*.fasta")])
        if file_path:
            try:
                result = self.result_text.get("1.0", tk.END)
                with open(file_path, "w") as file:
                    if file_path.endswith(".fasta"):
                        file.write(">Translated_Sequence\n")
                    file.write(result)
                messagebox.showinfo("Exported", "Result exported successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Could not export file: {e}")


if __name__ == "__main__":
    root = tk.Tk()
    app = TranslatorGUI(root)
    root.mainloop()
