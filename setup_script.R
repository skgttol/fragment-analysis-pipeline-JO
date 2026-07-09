# ==============================================================================
# MASTER PROJECT SETUP & ONBOARDING PIPELINE (USETHIS EDITION)
# ==============================================================================

# --- CONFIGURATION TOGGLES ---
force_update_scripts <- FALSE 

# Update with your specific GitHub details
github_owner <- "skgttol"
github_repo  <- "fragment-analysis-pipeline-JO"
github_repo_full <- paste0(github_owner, "/", github_repo) # "Owner/Repo"

# Define the direct zip source and target locations
zip_url   <- sprintf("https://github.com/%s/%s/archive/refs/heads/main.zip", github_owner, github_repo)
temp_zip  <- tempfile(fileext = ".zip")

# Ensure usethis is installed
if (!requireNamespace("usethis", quietly = TRUE)) {
  message("Installing required 'usethis' dependency...")
  install.packages("usethis", type = "binary")
}

# Ensure usethis is installed
if (!requireNamespace("rstudioapi", quietly = TRUE)) {
  message("Installing required 'rstudioapi' dependency...")
  install.packages("rstudioapi", type = "binary")
}

# ==============================================================================
# STEP 1: CONSOLIDATED PACKAGE INSTALLATION & CRAN BINARY FALLBACK
# ==============================================================================
message("\n--- PHASE 2: INITIALIZING PIPELINE PACKAGES ---")

# 1. Dynamically check R version to determine the best repository
current_r_version <- getRversion()

if (current_r_version < "4.4.0") {
  message(sprintf("Detected R version %s. Routing to Posit Package Manager for compatible binaries...", current_r_version))
  options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"))
} else {
  message(sprintf("Detected R version %s. Using default CRAN repository...", current_r_version))
  # Ensures a default repo is set so the script doesn't hang on a fresh install
  if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
}

# Comprehensive master package list across all scripts
required_pkgs <- c(
  "boot", "broom", "broom.mixed", "cowplot", "data.table", "datawizard", 
  "dendextend", "emmeans", "future", "GGally", "ggeffects", "ggh4x", 
  "ggnewscale", "ggplot2", "ggpubr", "ggrepel", "ggridges", "ggstats", 
  "grid", "gridExtra", "gtable", "here", "htmltools", "insight", 
  "janitor", "kableExtra", "knitr", "lme4", "lmerTest", "logr", "magick", "Matrix",
  "MASS", "openxlsx", "patchwork", "plotly", "ragg", "RColorBrewer", 
  "readxl", "reshape2", "rmarkdown", "scales", "segmented", "tidyverse", 
  "writexl", "yaml"
)

# Identify uninstalled packages
missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]

if (length(missing_pkgs) > 0) {
  # Dynamically check for any version of Rtools environment setups
  has_rtools <- any(grepl("RTOOLS", names(Sys.getenv()))) || nzchar(Sys.which("make"))
  
  if (has_rtools) {
    message("Rtools detected. Installing missing analysis packages...")
    install.packages(missing_pkgs)
  } else {
    message("Rtools NOT found. Installing pre-compiled binaries...")
    install.packages(missing_pkgs, type = getOption("pkgType")) 
  }
} else {
  message("✅ All required pipelines and dependencies are successfully installed.")
}

# Load packages cleanly into memory
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ============================================================================== 
# STEP 2: INTERACTIVE POPUPS FOR LOCATION AND PROJECT NAME
# ============================================================================== 
message("Waiting for user input...")

# 1. Popup to select the Parent Directory
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  parent_dir <- rstudioapi::selectDirectory(
    caption = "Select the parent directory where your project folder will be created"
  )
} else {
  parent_dir <- readline(prompt = "Enter the full path for the project parent directory: ")
}

if (is.null(parent_dir) || !nzchar(parent_dir)) {
  stop("Project setup cancelled: No directory was selected.", call. = FALSE)
}

parent_dir <- normalizePath(parent_dir, winslash = "/", mustWork = TRUE)

# 2. Popup to enter the Project Name (Using native RStudio prompt with console fallback)
project_name <- ""

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  # Native RStudio modal prompt (guaranteed top-level focus)
  project_name <- rstudioapi::showPrompt(
    title   = "Project Name",
    message = "Enter a name for your new project folder:",
    default = "My_Analysis_Pipeline"
  )
} else {
  # Fallback if running from a standard terminal/GUI outside RStudio
  project_name <- readline(prompt = "Enter project name [default: My_Analysis_Pipeline]: ")
  if (!nzchar(project_name)) project_name <- "My_Analysis_Pipeline"
}

# Halt if user clicks 'Cancel' on the prompt
if (is.null(project_name)) {
  stop("Project setup cancelled: No project name entered.", call. = FALSE)
}

# Format folder name safely (strip extra spaces and replace with underscores)
project_name <- gsub("\\s+", "_", trimws(project_name))
project_dir  <- file.path(parent_dir, project_name)

# ==============================================================================
# STEP 3 & 4: DOWNLOAD FULL REPOSITORY ZIP (COLLISION PROOF)
# ==============================================================================
message("\n--- PHASE 1: DOWNLOADING WORKSPACE SNAPSHOT ---")

# We look for a 'scripts/' folder inside the new path to see if it's already set up
scripts_dir <- file.path(project_dir, "analysis-scripts")
is_empty    <- !dir.exists(scripts_dir) || length(list.files(scripts_dir)) == 0

if (is_empty || force_update_scripts) {
  message("--> Requesting secure archive download from GitHub...")
  
  # Ensure target project root exists before unzipping
  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    # Use your proven auth headers to fetch the full archive directly
    dl_status <- download.file(
      url = zip_url, 
      destfile = temp_zip, 
      mode = "wb",          # Prevents zip file corruption on Windows
      quiet = TRUE
    )
    
    if (dl_status != 0) stop("Download rejected or connection timed out.")
    
    # Extract the archive into the newly generated folder
    unzip(zipfile = temp_zip, exdir = project_dir)
    
    # GitHub names branch archives as "RepoName-BranchName"
    expected_wrapper <- file.path(project_dir, paste0(github_repo, "-main"))
    
    if (dir.exists(expected_wrapper)) {
      # Lift all files and nested subdirectories out of the wrapper folder
      top_level_items <- list.files(expected_wrapper, full.names = TRUE, all.files = TRUE, no.. = TRUE)
      
      # Recursively move everything up directly into the project directory
      file.copy(from = top_level_items, to = project_dir, overwrite = TRUE, recursive = TRUE)
      
      # Purge temporary wrapper directory and download artifacts
      unlink(expected_wrapper, recursive = TRUE)
      unlink(temp_zip)

      message("✅ Alls cript trees and templates deployed successfully.")
    } else {
      stop("Extraction structure mismatch. Expected branch folder not found.")
    }
    
  }, error = function(e) {
    unlink(temp_zip, showWarnings = FALSE)
    stop(paste("\n❌ ERROR: Secure workspace download failed.", e$message), call. = FALSE)
  })
  
} else {
  message("--> Project already populated. Skipping download. (Set force_update_scripts <- TRUE to update)")
}


# ==============================================================================
# STEP 5: SETUP AUTO-LAUNCH & OPEN NEW PROJECT INSTANCE
# ==============================================================================

# Specify the relative path of the file you want to automatically open inside the new project
target_script_to_open <- "master_script.R"

# 1. Write a self-deleting .Rprofile into the new project root with bulletproof timing
rprofile_path <- file.path(project_dir, ".Rprofile")

rprofile_content <- c(
  "# Temporary startup automation generated by setup pipeline",
  "if (interactive()) {",
  "  # Use RStudio's session initialization hook",
  "  setHook('rstudio.sessionInit', function(newSession) {",
  "    # Use a tiny 0.5 second delay to ensure graphical panes finish rendering",
  "    if (requireNamespace('rstudioapi', quietly = TRUE)) {",
  "      # Force absolute path construction within the newly established workspace",
  sprintf("      target_file <- file.path(getwd(), '%s')", target_script_to_open),
  "      if (file.exists(target_file)) {",
  "        # Instruct RStudio IDE directly to navigate to and open the document",
  "        rstudioapi::navigateToFile(target_file)",
  "        message('--> Master pipeline script auto-launched.')",
  "      } else {",
  "        warning('Auto-launch failed: Target script could not be located at ', target_file)",
  "      }",
  "    }",
  "    # Self-destruct: remove this .Rprofile so it only runs once on initial setup",
  "    unlink('.Rprofile', force = TRUE)",
  "  }, action = 'append')",
  "}"
)

writeLines(rprofile_content, con = rprofile_path)

# 2. Write standard RStudio project configuration file
rproj_path <- file.path(project_dir, paste0(project_name, ".Rproj"))
if (!file.exists(rproj_path)) {
  rproj_config <- c(
    "Version: 1.0", "", "RestoreWorkspace: No", "SaveWorkspace: No",
    "AlwaysSaveHistory: Default", "", "EnableCodeIndexing: Yes",
    "UseSpacesForTab: Yes", "NumSpacesForTab: 2", "Encoding: UTF-8"
  )
  writeLines(rproj_config, con = rproj_path)
}

# 3. Open the new project folder in a fresh RStudio window instance
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  rstudioapi::openProject(path = project_dir, newSession = TRUE)
}

# 4. Cleanly halt automated setup execution in the original window
message(
  paste0(
    "\n=====================================================================\n",
    "SETUP COMPLETE & HALTED\n",
    "=====================================================================\n",
    "Your project workspace has been successfully built at:\n", project_dir, "\n\n",
    "A new RStudio window is launching with your master script opened!\n",
    "====================================================================="
  ),
  call. = FALSE
)
