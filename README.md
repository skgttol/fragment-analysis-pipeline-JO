***

# 🧬 Advanced Fragment Analysis Pipeline

Welcome to the Advanced Fragment Analysis Pipeline. This repository contains an automated, end-to-end R-based workflow for processing, modeling, and visualizing fragment analysis data (e.g., CAG repeat sizing, somatic instability). 

This pipeline is designed to be highly reproducible. You only need to download a single **Master Run Script** and your templates; the script will handle syncing the rest of the code directly from this repository.

---

## 📋 Prerequisites
Before you begin, ensure you have the following installed on your computer:
1. **[R](https://cran.r-project.org/)** (v4.1.0 or higher)
2. **[RStudio](https://posit.co/download/rstudio-desktop/)**

---

## 🚀 Step 1: Initial Download & Folder Setup

To keep your analysis organized and reproducible, we use a strict folder structure and an RStudio Project (`.Rproj`).

### 1. Download the Core Files
Download the following three files from this GitHub repository and place them into an empty folder on your computer (e.g., `Documents/Fragment_Analysis_Exp1/`):
* `run_pipeline.R` (The Master Control Script)
* `config_template.yml` (Rename this to `config.yml` once downloaded)
* `platemap_template.xlsx` (Rename as needed for your experiment)

### 2. Create the Folder Structure
Inside your new project folder, create the following sub-folders. *(Note: The pipeline will look for these specific folders!)*

```text
📁 Your_Project_Folder/
│
├── 📄 run_pipeline.R          <- You downloaded this
├── 📄 config.yml              <- You downloaded & renamed this
├── 📄 platemap.xlsx           <- You downloaded & renamed this
│
├── 📁 1_raw_fsa/              <- Drop your raw .fsa files here
├── 📁 2_metadata/            <- Move your platemap.xlsx here
└── 📁 3_romeo_output/      <- (Empty for now) External program output goes here
```

### 3. Initialize the RStudio Project
1. Open RStudio.
2. Go to **File** > **New Project** > **Existing Directory**.
3. Select your `Your_Project_Folder` and click **Create Project**.
*(This creates a `.Rproj` file. Always double-click this file to open your project in the future. It ensures R always knows exactly where your files are).*

---

## ⚙️ Step 2: Configuration & Platemaps

1. **Fill out the Platemap:** Open your platemap Excel file in `2_metadata/` and ensure your samples are mapped to their correct wells and plates.
2. **Edit `config.yml`:** Open `config.yml` in RStudio or a text editor. Update the parameters to match your experiment:
   * Point the paths to your `1_raw_fsa/` folder and your platemap file.
   * Define your primary and secondary grouping variables (e.g., Genotype, Clone).
   * Define your analysis windows and thresholds.

---

## 🏃‍♂️ Step 3: Running the Pipeline (Phase 1)

Before the main R pipeline can run, we need to generate settings for Romeo.

1. Open `run_pipeline.R` in RStudio.
2. At the top of the script in the **Control Panel**, set the toggles as follows:
   ```r
   RUN_00_SETTINGS    <- TRUE   # Turn this ON for the first run
   RUN_01_PROCESS     <- FALSE  
   RUN_02_PLOTS       <- FALSE  
   # ... (set the rest to FALSE)
   ```
3. Run the entire `run_pipeline.R` script.
   * *The script will securely connect to GitHub, download the latest analysis scripts in the background, and generate your custom settings file.*
4. Open your external sizing program, import the generated settings file, and run your raw `.txt` files from Genemapper through it.
5. Export the sizing results (Excel file) from the external program and save it into the **`3_romeo_output/`** folder.

---

## 📊 Step 4: Running the Pipeline (Phase 2)

Now that your data has been externally sized, you are ready to run the full automated analysis.

1. Update your `config.yml` so that `external_data_path` points to the new Excel file in your `3_romeo_output/` folder.
2. Open `run_pipeline.R` in RStudio.
3. Update the **Control Panel** toggles:
   ```r
   RUN_00_SETTINGS    <- FALSE  # Turn this OFF now
   
   RUN_01_PROCESS     <- TRUE   # Turn all core scripts ON
   RUN_02_PLOTS       <- TRUE  
   RUN_03_CORRELATION <- TRUE  
   RUN_04_FIGURES     <- TRUE  
   RUN_05_TRACES      <- TRUE  
   ```
4. Run the entire `run_pipeline.R` script.

---

## 📂 Understanding Your Output

The pipeline will automatically create an `Output/` directory stamped with today's date and time. Inside, you will find:

* **`01_QC_plots/`**: Diagnostic plots checking for biological and technical outliers.
* **`02_descriptive_plots/`**: Group averages, individual replicate trajectories, and heatmaps.
* **`03_model_plots/`**: Linear mixed-effects model outputs, confidence intervals, and waterfall plots.
* **`04_correlation/`**: PCA biplots, hierarchical clustering dendrograms, and correlation matrices.
* **`combined_figures/`**: Publication-ready, multi-panel stitched figures.
* **`trace_plots/`**: Multi-page PDF trace reports showing start vs. end expansions, overlaid with your baseline and thresholds.
* **Data Tables**: Cleaned, wide-format Excel sheets ready for GraphPad Prism or custom querying.

*(Note: Future versions of this pipeline will compile these outputs into a single, interactive HTML dashboard for easy browsing).*

---
### 💡 Troubleshooting
* **Error: "File not found"**: Ensure you have opened RStudio by double-clicking the `.Rproj` file, and verify your folder names exactly match what is written in `config.yml`.
* **Script fails to download**: Check your internet connection. If the repository is private, ensure your Personal Access Token (PAT) in `run_pipeline.R` is valid and hasn't expired.
