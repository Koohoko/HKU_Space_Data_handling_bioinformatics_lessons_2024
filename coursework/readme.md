# Scientific Skills Bioinformatics Exercise (30%)

**Module Title**: Scientific Skills  
**Program**: MSc Biomedical Science  
**Assessment Weighting**: 30% of Module Marks  

---

## Exercise Overview
This exercise introduces students to bioinformatics data handling, focusing on data cleaning, sequence analysis, and reproducibility. Students will work in a pre-configured Google Colab notebook with guided sections to facilitate learning without requiring prior programming experience.

---

## Learning Objectives
1. Perform basic bioinformatics data cleaning and analysis on sequence data.
2. Calculate GC content and codon usage within specific regions of SARS-CoV-2 sequences.
3. Apply reproducibility and good organization practices in bioinformatics workflows.

---

## Instructions for Accessing the Colab Notebook
1. **Access the Notebook**: Open the link provided on Moodle to access the pre-configured Google Colab notebook.
2. **Notebook Overview**: The notebook is divided into sections with instructions and explanations for each step.
3. **Running Cells**: Click on each code cell and press the “Run” button to execute the code and view results.
4. **Follow Along**: Each code block is accompanied by an explanation, so follow along and read the comments carefully.

---

## Exercise Outline and Assessment Tasks

### 1. Introduction to Bioinformatics Data Skills
- **Background**: Brief overview of bioinformatics data formats, such as FASTA, and the importance of reproducibility. Students will work with SARS-CoV-2 sequences, focusing on calculating coverage and identifying specific genomic regions.

### 2. Data Cleaning and Quality Control
- **Tasks**:
  1. Load the provided SARS-CoV-2 FASTA file (20 sequences, each 29,903 bp) into the Colab notebook.
  2. Filter out sequences with coverage below 85% (counting only A, C, T, and G bases).
  3. Summarize the cleaning process and explain how data quality impacts analysis.
- **Code Block**: The Colab notebook includes a code block that loads and filters the sequences based on coverage. Students need to run the cell and observe the output.

### 3. Sequence Analysis
- **Tasks**:
  1. Calculate GC content for each of two randomly selected sequences.
  2. Extract the spike gene region (positions 21,563 to 25,384) for both sequences.
  3. Calculate the codon usage for one of the extracted sequences.
- **Guided Code Blocks**: Each analysis task has its own code block with comments explaining what each line does. Students simply run the cells and observe the results.

### 4. Reproducibility and Data Organization
- **Tasks**:
  1. Follow best practices for naming and organizing files within the notebook.
  2. Ensure reproducibility by setting a random seed for the code block that selects sequences.
  3. Observe comments within the notebook, explaining how reproducibility is maintained.
- **Guided Example**: Students will be guided to set up a structured notebook document in Colab and see an example of setting a random seed to make results reproducible.

### 5. Reflection on Bioinformatics Data Skills
- Reflect on the importance of reproducibility and good data organization practices in bioinformatics workflows.

---

## Submission Guidelines
The report should include:

1. **Introduction (10 marks)**
   - Overview of bioinformatics data skills and significance.
2. **Data Cleaning and Quality Control Summary (20 marks)**
   - Explanation of the data cleaning process based on coverage and its importance.
3. **Sequence Analysis (35 marks)**
   - Results and interpretation of GC content and codon usage analysis (Are the results the same for the two sequences? Interpret your observation).
4. **Notebook Organization and Reproducibility (20 marks)**
   - A print of the organized notebook which you used for generating the results, with readable code and ensuring the results are reproducible.
5. **Reflection (15 marks)**
   - Reflection on the importance of reproducibility and data organization.

---