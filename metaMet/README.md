# metaMet

## Data Preprocessing Pipeline

This repository contains the data preprocessing pipeline for metaMet. Below are the prerequisites and installation instructions to set up the pipeline.

---

## Prerequisites

### BRENDA JSON File

1. **Download the BRENDA JSON File:**
   - Visit the official BRENDA download page: [BRENDA Download](https://www.brenda-enzymes.org/download.php).
   - Download the JSON file.

2. **Place the JSON File:**
   - Move the downloaded JSON file to the following directory:
     ```
     metaMet/data/raw/brenda/
     ```
   - **Important:** Ensure there is only one JSON file in the `brenda` folder.

---

### BioCyc Databases (MetaCyc)

1. **Download the BioCyc Flatfile Databases:**
   - Visit the BioCyc download page: [BioCyc Flatfiles](https://bioinformatics.ai.sri.com/ecocyc/dist/flatfiles-52983746/).
   - Download the Tier 1 and Tier 2 databases in flatfile format.
   - **Note:** You will need to contact BioCyc for credentials to access the files.

2. **Place the Databases:**
   - Move the downloaded directory to:
     ```
     metaMet/data/raw/metacyc/
     ```
   - **Important:** Ensure there is only one directory in the `metacyc` folder, and its name **must not** contain a period (`.`).

---

## Installation

To set up the pipeline, follow these steps:

1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/metaMet.git