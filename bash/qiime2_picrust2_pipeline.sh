#!/usr/bin/env bash
#===============================================================================
# Script: qiime2_picrust2_pipeline.sh
# Author: Tiago Mota
# Date: 2025-09-08
# Description:
#   QIIME2 16S rRNA processing pipeline with optional PICRUSt2 functional
#   prediction. Designed to run standalone or within SLURM jobs.
#
# Usage:
#   ./qiime2_picrust2_pipeline.sh [OPTIONS]
#
# SILVA Database Options:
#   - If primers are provided (--forward-primer and --reverse-primer):
#     The database will be trimmed to the specified region and a region-specific
#     classifier will be built (classifier.qza). This is faster for classification.
#   - If no primers are provided:
#     The full-length database will be used (classifier_all.qza). This is more
#     comprehensive but slower. A message will be printed to inform the user.
#===============================================================================

set -euo pipefail

#-----------------------------
# Global Configuration
#-----------------------------
SCRIPT_NAME="$(basename "$0")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default values
DENOISE="deblur"
PICRUST="false"
THREADS=4
QIIME_ENV="qiime2-amplicon-2025.7"
PICRUST_ENV="picrust2"
SILVA_VERSION="138.2"
SILVA_TARGET="SSURef_NR99"
FORWARD_PRIMER=""  # Empty by default
REVERSE_PRIMER=""  # Empty by default
SILVA_DIR=""       # Will be set based on work-dir

#-----------------------------
# Functions
#-----------------------------
usage() {
    cat << EOF
Usage: ${SCRIPT_NAME} [OPTIONS]

Required:
  -w, --work-dir <path>          Working directory (absolute path)
  -s, --seqs-dir <path>          Directory containing FASTQ files
  -f, --forward-pattern <str>    Pattern for forward reads (e.g., _R1, _1.fq.gz)
  -r, --reverse-pattern <str>    Pattern for reverse reads (e.g., _R2, _2.fq.gz)
  -n, --n-char <int>             Characters to trim for sample ID (excludes patterns)
  -m, --o-manifest <path>        Output manifest file path (without extension)
  -q, --o-qiime <path>           Output directory for QIIME2 results

Optional:
  -d, --d-denoise <method>       Denoising method: deblur (default) or dada2
  -p, --p-picrust <bool>         Run PICRUSt2 functional prediction (true/false)
  -t, --threads <int>            Number of threads (default: 4)
  --qiime-env <name>             Conda environment name for QIIME2 (default: qiime2-amplicon-2025.7)
  --picrust-env <name>           Conda environment name for PICRUSt2 (default: picrust2)
  --silva-dir <path>             SILVA database directory (default: work-dir/../SILVA)
  --silva-version <ver>          SILVA version (default: 138.2)
  
  Primer Options (for region-specific classification):
  --forward-primer <seq>         Forward primer sequence (e.g., CCTACGGGNGGCWGCAG)
  --reverse-primer <seq>         Reverse primer sequence (e.g., GGACTACNVGGGTWTCTAAT)
                                 IMPORTANT: Both primers must be provided together.
                                 If provided, a region-specific classifier will be built.
                                 If omitted, the full-length classifier will be used.
  
  -h, --help                     Show this help message

Examples:
  # Full-length classifier (no primers provided)
  ${SCRIPT_NAME} \\
    --work-dir /project/microbiome \\
    --seqs-dir /project/data/reads \\
    --forward-pattern _R1 \\
    --reverse-pattern _R2 \\
    --n-char 5 \\
    --o-manifest results/manifest \\
    --o-qiime results/qiime2

  # Region-specific classifier (V3-V4 with primers)
  ${SCRIPT_NAME} \\
    --work-dir /project/microbiome \\
    --seqs-dir /project/data/reads \\
    --forward-pattern _R1 \\
    --reverse-pattern _R2 \\
    --n-char 5 \\
    --o-manifest results/manifest \\
    --o-qiime results/qiime2 \\
    --forward-primer CCTACGGGNGGCWGCAG \\
    --reverse-primer GGACTACNVGGGTWTCTAAT

EOF
}

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

error() {
    echo "[ERROR] $*" >&2
    exit 1
}

check_command() {
    if ! command -v "$1" &> /dev/null; then
        error "Command '$1' not found. Please install or load the appropriate module."
    fi
}

check_conda_env() {
    local env_name="$1"
    if ! conda env list | grep -q "^${env_name} "; then
        error "Conda environment '$env_name' not found. Available environments:"
        conda env list
        return 1
    fi
    log "Conda environment '$env_name' found"
}

activate_env() {
    local env_name="$1"
    log "Activating conda environment: $env_name"
    
    # Initialize conda for bash
    if [[ -f "$CONDA_PREFIX/etc/profile.d/conda.sh" ]]; then
        source "$CONDA_PREFIX/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    else
        error "Cannot find conda.sh. Please ensure conda is installed and accessible."
    fi
    
    conda activate "$env_name" || error "Failed to activate environment: $env_name"
    log "Activated environment: $env_name (Python: $(python --version 2>&1))"
}

validate_args() {
    local ok=true

    # Check directories
    if [[ -z "$WDIR" ]]; then
        error "--work-dir is required"
    fi
    if [[ ! -d "$WDIR" ]]; then
        error "Directory '$WDIR' does not exist"
    fi

    if [[ -z "$SDIR" ]]; then
        error "--seqs-dir is required"
    fi
    if [[ ! -d "$SDIR" ]]; then
        error "Directory '$SDIR' does not exist"
    fi

    # Check patterns
    if [[ -z "$FORWARD_PATTERN" ]]; then
        error "--forward-pattern is required and cannot be empty"
    fi
    if [[ -z "$REVERSE_PATTERN" ]]; then
        error "--reverse-pattern is required and cannot be empty"
    fi

    # Check N_CHAR
    if [[ -z "$N_CHAR" ]]; then
        error "--n-char is required"
    fi
    if ! [[ "$N_CHAR" =~ ^[0-9]+$ ]]; then
        error "--n-char must be a non-negative integer (got '$N_CHAR')"
    fi

    # Check manifest output
    if [[ -z "$OUT_MANI" ]]; then
        error "--o-manifest is required and cannot be empty"
    else
        local outmdir=$(dirname "$OUT_MANI")
        mkdir -p "$outmdir" || error "Cannot create output directory '$outmdir'"
    fi

    # Check QIIME output
    if [[ -z "$OUT_QIIME" ]]; then
        error "--o-qiime is required and cannot be empty"
    else
        mkdir -p "$OUT_QIIME" || error "Cannot create QIIME output directory '$OUT_QIIME'"
    fi

    # Validate denoise method
    if [[ "$DENOISE" != "deblur" && "$DENOISE" != "dada2" ]]; then
        log "Warning: --d-denoise should be 'deblur' or 'dada2'. Using default: 'deblur'."
        DENOISE="deblur"
    fi

    # Validate PICRUST
    if [[ "$PICRUST" =~ ^(true|True|T|TRUE|false|False|F|FALSE)$ ]]; then
        if [[ "$PICRUST" =~ ^(true|True|T|TRUE)$ ]]; then
            PICRUST="true"
        else
            PICRUST="false"
        fi
    else
        log "Warning: --p-picrust should be true or false. Using default: 'false'."
        PICRUST="false"
    fi

    # Validate threads
    if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
        error "--threads must be a positive integer (got '$THREADS')"
    fi

    # Validate primer logic
    if [[ -n "$FORWARD_PRIMER" && -z "$REVERSE_PRIMER" ]]; then
        error "Both --forward-primer and --reverse-primer must be provided together. If you want to use the full-length classifier, omit both."
    fi
    if [[ -z "$FORWARD_PRIMER" && -n "$REVERSE_PRIMER" ]]; then
        error "Both --forward-primer and --reverse-primer must be provided together. If you want to use the full-length classifier, omit both."
    fi

    # Set SILVA_DIR if not provided
    if [[ -z "$SILVA_DIR" ]]; then
        SILVA_DIR="$WDIR/../SILVA"
    fi
    mkdir -p "$SILVA_DIR"

    return 0
}

setup_silva_database() {
    log "=== Setting up SILVA database ==="
    
    # Determine which classifier to use based on primer presence
    local classifier="${SILVA_DIR}/classifier.qza"
    local classifier_all="${SILVA_DIR}/classifier_all.qza"
    
    # Check if primers were provided
    if [[ -n "$FORWARD_PRIMER" && -n "$REVERSE_PRIMER" ]]; then
        log "🔬 Primers provided: Forward='$FORWARD_PRIMER', Reverse='$REVERSE_PRIMER'"
        log "Building region-specific classifier (this is faster for classification)"
        
        if [[ -f "$classifier" ]]; then
            log "Region-specific classifier already exists: $classifier"
            return 0
        fi
        
        build_region_specific_classifier
    else
        log "⚠️  No primers provided. Using full-length SILVA database."
        log "ℹ️  Since no specific region was selected, all variable regions will be used."
        log "ℹ️  Classification will be comprehensive but may take longer."
        log "ℹ️  To build a region-specific classifier (faster), provide --forward-primer and --reverse-primer."
        
        if [[ -f "$classifier_all" ]]; then
            log "Full-length classifier already exists: $classifier_all"
            return 0
        fi
        
        build_full_length_classifier
    fi
}

build_region_specific_classifier() {
    local original_dir=$(pwd)
    cd "$SILVA_DIR" || error "Cannot cd to $SILVA_DIR"
    
    log "=== Building region-specific SILVA classifier ==="
    log "Region: Forward='$FORWARD_PRIMER', Reverse='$REVERSE_PRIMER'"
    
    # Get Silva data if not present
    if [[ ! -f "silva-${SILVA_VERSION}-ssu-nr99-seqs.qza" ]]; then
        log "Downloading SILVA database (version ${SILVA_VERSION})..."
        qiime rescript get-silva-data \
            --p-version "${SILVA_VERSION}" \
            --p-target "${SILVA_TARGET}" \
            --p-include-species-labels \
            --o-silva-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs.qza" \
            --o-silva-taxonomy "silva-${SILVA_VERSION}-ssu-nr99-tax.qza" \
            || error "Failed to download SILVA data"
    fi
    
    # Clean sequences
    if [[ ! -f "silva-${SILVA_VERSION}-ssu-nr99-seqs-cleaned.qza" ]]; then
        log "Cleaning sequences..."
        qiime rescript cull-seqs \
            --i-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs.qza" \
            --o-clean-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs-cleaned.qza" \
            || error "Failed to clean sequences"
    fi
    
    # Filter unwanted taxa
    if [[ ! -f "silva-${SILVA_VERSION}-ssu-nr99-seqs-filt.qza" ]]; then
        log "Filtering unwanted taxa..."
        qiime rescript filter-seqs-length-by-taxon \
            --i-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs-cleaned.qza" \
            --i-taxonomy "silva-${SILVA_VERSION}-ssu-nr99-tax.qza" \
            --p-labels Archaea Bacteria Eukaryota Mitochondria Chloroplast Unknown NA metagenome unidentified \
            --p-min-lens 900 1200 5000 5000 5000 5000 5000 5000 5000 \
            --o-filtered-seqs "silva-${SILVA_VERSION}-ssu-nr99-seqs-filt.qza" \
            --o-discarded-seqs "silva-${SILVA_VERSION}-ssu-nr99-seqs-discard.qza" \
            || error "Failed to filter sequences"
    fi
    
    # Dereplicate
    if [[ ! -f "silva-${SILVA_VERSION}-ssu-nr99-tax-derep-super.qza" ]]; then
        log "Dereplicating sequences..."
        qiime rescript dereplicate \
            --i-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs-filt.qza" \
            --i-taxa "silva-${SILVA_VERSION}-ssu-nr99-tax.qza" \
            --p-mode 'uniq' \
            --p-threads "$THREADS" \
            --o-dereplicated-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs-derep-super.qza" \
            --o-dereplicated-taxa "silva-${SILVA_VERSION}-ssu-nr99-tax-derep-super.qza" \
            || error "Failed to dereplicate sequences"
    fi
    
    # Extract reads with the provided primers
    log "Extracting reads for region defined by primers..."
    qiime feature-classifier extract-reads \
        --i-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs-filt.qza" \
        --p-f-primer "$FORWARD_PRIMER" \
        --p-r-primer "$REVERSE_PRIMER" \
        --p-trunc-len 120 \
        --p-min-length 100 \
        --p-max-length 400 \
        --o-reads "ref-seqs_region-specific.qza" \
        || error "Failed to extract reads with provided primers"
    
    # Build region-specific classifier
    log "Building region-specific classifier..."
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads "ref-seqs_region-specific.qza" \
        --i-reference-taxonomy "silva-${SILVA_VERSION}-ssu-nr99-tax-derep-super.qza" \
        --o-classifier "classifier.qza" \
        || error "Failed to build region-specific classifier"
    
    cd "$original_dir" || error "Cannot return to $original_dir"
    log "✅ Region-specific classifier built successfully: $SILVA_DIR/classifier.qza"
}

build_full_length_classifier() {
    local original_dir=$(pwd)
    cd "$SILVA_DIR" || error "Cannot cd to $SILVA_DIR"
    
    log "=== Building full-length SILVA classifier ==="
    log "ℹ️  This will classify against the entire SILVA database."
    log "ℹ️  Expect longer processing time compared to region-specific classification."
    
    # Get Silva data if not present
    if [[ ! -f "silva-${SILVA_VERSION}-ssu-nr99-seqs.qza" ]]; then
        log "Downloading SILVA database (version ${SILVA_VERSION})..."
        qiime rescript get-silva-data \
            --p-version "${SILVA_VERSION}" \
            --p-target "${SILVA_TARGET}" \
            --p-include-species-labels \
            --o-silva-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs.qza" \
            --o-silva-taxonomy "silva-${SILVA_VERSION}-ssu-nr99-tax.qza" \
            || error "Failed to download SILVA data"
    fi
    
    # Clean sequences
    if [[ ! -f "silva-${SILVA_VERSION}-ssu-nr99-seqs-cleaned.qza" ]]; then
        log "Cleaning sequences..."
        qiime rescript cull-seqs \
            --i-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs.qza" \
            --o-clean-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs-cleaned.qza" \
            || error "Failed to clean sequences"
    fi
    
    # Build full-length classifier
    log "Building full-length classifier (this may take a while)..."
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads "silva-${SILVA_VERSION}-ssu-nr99-seqs-cleaned.qza" \
        --i-reference-taxonomy "silva-${SILVA_VERSION}-ssu-nr99-tax.qza" \
        --o-classifier "classifier_all.qza" \
        || error "Failed to build full-length classifier"
    
    cd "$original_dir" || error "Cannot return to $original_dir"
    log "✅ Full-length classifier built successfully: $SILVA_DIR/classifier_all.qza"
}

create_manifest() {
    log "=== Creating manifest file ==="
    cd "$WDIR" || error "Cannot cd to $WDIR"
    
    # Get file lists
    local file_list=($SDIR/*)
    local for_seqs=($(printf '%s\n' "${file_list[@]}" | grep "$FORWARD_PATTERN" || true))
    local rev_seqs=($(printf '%s\n' "${file_list[@]}" | grep "$REVERSE_PATTERN" || true))
    
    if [[ ${#for_seqs[@]} -eq 0 ]] || [[ ${#rev_seqs[@]} -eq 0 ]]; then
        error "No files found matching patterns: forward='$FORWARD_PATTERN', reverse='$REVERSE_PATTERN'"
    fi
    
    if [[ ${#for_seqs[@]} -ne ${#rev_seqs[@]} ]]; then
        error "Number of forward (${#for_seqs[@]}) and reverse (${#rev_seqs[@]}) files do not match."
    fi
    
    # Extract replicate IDs
    local rep_id=()
    for f in "${for_seqs[@]}"; do
        local id=$(basename "$f" | sed -e "s/${FORWARD_PATTERN}//")
        rep_id+=("$id")
    done
    
    # Write manifest
    {
        echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"
        for ((i=0; i<${#rep_id[@]}; i++)); do
            echo -e "${rep_id[i]}\t${for_seqs[i]}\t${rev_seqs[i]}"
        done
    } > "${OUT_MANI}.tsv"
    
    log "Manifest created: ${OUT_MANI}.tsv (${#rep_id[@]} samples)"
}

extract_from_qza() {
    local qza_file="$1"
    local file_pattern="$2"
    
    if [[ ! -f "$qza_file" ]]; then
        error "QZA file not found: $qza_file"
    fi
    
    local folder_name=$(unzip -Z -1 "$qza_file" 2>/dev/null | head -1 | sed 's/\/[^/]*$//')
    if [[ -z "$folder_name" ]]; then
        error "Cannot determine folder structure in $qza_file"
    fi
    
    local extracted_file="$OUT_QIIME/${folder_name}/data/${file_pattern}"
    if [[ ! -f "$extracted_file" ]]; then
        log "Extracting $qza_file..."
        unzip -q "$qza_file" -d "$OUT_QIIME/" || error "Failed to extract $qza_file"
    fi
    
    echo "$extracted_file"
}

get_optimal_trunc_length() {
    local qzv_file="$1"
    local folder_name=$(unzip -Z -1 "$qzv_file" 2>/dev/null | head -1 | sed 's/\/[^/]*$//')
    
    if [[ -z "$folder_name" ]]; then
        error "Cannot determine folder structure in $qzv_file"
    fi
    
    local tsv_file="$OUT_QIIME/${folder_name}/data/forward-seven-number-summaries.tsv"
    if [[ ! -f "$tsv_file" ]]; then
        unzip -q "$qzv_file" -d "$OUT_QIIME/" || error "Failed to extract $qzv_file"
    fi
    
    local target_col=$(awk -F'\t' 'NR==2 {for(i=2;i<=NF;i++) if($i+0 < 9000) {print i-2; exit}}' "$tsv_file")
    target_col=${target_col:-0}
    
    echo "$target_col"
}

run_qiime2() {
    log "=== Starting QIIME2 processing with environment: $QIIME_ENV ==="
    
    # Activate QIIME2 environment
    activate_env "$QIIME_ENV"
    
    # Import sequences
    local demux="$OUT_QIIME/paired-end-demux.qza"
    if [[ ! -f "$demux" ]]; then
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path "${OUT_MANI}.tsv" \
            --output-path "$demux" \
            --input-format PairedEndFastqManifestPhred33V2 \
            || error "Failed to import sequences"
    fi
    
    # Denoising
    if [[ "$DENOISE" == "deblur" ]]; then
        run_deblur "$demux"
    else
        run_dada2 "$demux"
    fi
    
    # Taxonomy classification - determine which classifier to use
    log "=== Classifying taxonomy ==="
    local classifier=""
    if [[ -n "$FORWARD_PRIMER" && -n "$REVERSE_PRIMER" ]]; then
        classifier="$SILVA_DIR/classifier.qza"
        log "Using region-specific classifier (primers: $FORWARD_PRIMER / $REVERSE_PRIMER)"
    else
        classifier="$SILVA_DIR/classifier_all.qza"
        log "Using full-length classifier (no region specified)"
    fi
    
    if [[ ! -f "$classifier" ]]; then
        error "Classifier not found: $classifier"
    fi
    
    qiime feature-classifier classify-sklearn \
        --i-classifier "$classifier" \
        --i-reads "$OUT_DENOISE" \
        --p-confidence 0.97 \
        --p-n-jobs "$THREADS" \
        --o-classification "$OUT_TAXONOMY" \
        || error "Taxonomy classification failed"
    
    # Phylogenetic tree
    log "=== Building phylogenetic tree ==="
    if [[ ! -f "$OUT_ALIGNED" ]]; then
        qiime alignment mafft \
            --i-sequences "$OUT_DENOISE" \
            --o-alignment "$OUT_ALIGNED" \
            --p-n-threads "$THREADS" \
            || error "MAFFT alignment failed"
    fi
    
    if [[ ! -f "$OUT_MASKED" ]]; then
        qiime alignment mask \
            --i-alignment "$OUT_ALIGNED" \
            --o-masked-alignment "$OUT_MASKED" \
            || error "Masking failed"
    fi
    
    if [[ ! -f "$OUT_TREE" ]]; then
        qiime phylogeny fasttree \
            --i-alignment "$OUT_MASKED" \
            --o-tree "$OUT_TREE" \
            --p-n-threads "$THREADS" \
            || error "FastTree failed"
    fi
    
    # Deactivate QIIME2 environment
    conda deactivate
    
    log "QIIME2 processing complete!"
}

run_deblur() {
    local demux="$1"
    log "=== Running Deblur denoising ==="
    
    # Merge paired ends
    local joined="$OUT_QIIME/demux-joined.qza"
    if [[ ! -f "$joined" ]]; then
        qiime vsearch merge-pairs \
            --i-demultiplexed-seqs "$demux" \
            --o-merged-sequences "$joined" \
            --o-unmerged-sequences "$OUT_QIIME/demux-unmerged.qza" \
            || error "VSEARCH merge failed"
    fi
    
    # Quality filter
    local filtered="$OUT_QIIME/demux-joined-filtered.qza"
    if [[ ! -f "$filtered" ]]; then
        qiime quality-filter q-score \
            --i-demux "$joined" \
            --p-min-quality 30 \
            --o-filtered-sequences "$filtered" \
            --o-filter-stats "$OUT_QIIME/demux-joined-filter-stats.qza" \
            || error "Quality filtering failed"
    fi
    
    # Visualize
    local qzv="$OUT_QIIME/demux-joined-filtered.qzv"
    if [[ ! -f "$qzv" ]]; then
        qiime demux summarize \
            --i-data "$filtered" \
            --o-visualization "$qzv" \
            || error "Demux summarization failed"
    fi
    
    # Get optimal truncation length
    local trim_length=$(get_optimal_trunc_length "$qzv")
    log "Optimal truncation length: $trim_length"
    
    # Set output variables
    OUT_DENOISE="$OUT_QIIME/rep-seqs-deblur.qza"
    OUT_TABLE="$OUT_QIIME/table-deblur.qza"
    OUT_TAXONOMY="$OUT_QIIME/taxonomy-deblur.qza"
    OUT_ALIGNED="$OUT_QIIME/aligned-rep-seqs-deblur.qza"
    OUT_MASKED="$OUT_QIIME/masked-aligned-rep-seqs-deblur.qza"
    OUT_TREE="$OUT_QIIME/fasttree-tree-deblur.qza"
    
    # Run deblur
    if [[ ! -f "$OUT_DENOISE" ]]; then
        qiime deblur denoise-16S \
            --i-demultiplexed-seqs "$filtered" \
            --p-trim-length "$trim_length" \
            --p-sample-stats \
            --p-jobs-to-start "$THREADS" \
            --o-representative-sequences "$OUT_DENOISE" \
            --o-table "$OUT_TABLE" \
            --o-stats "$OUT_QIIME/deblur-stats.qza" \
            || error "Deblur failed"
    fi
}

run_dada2() {
    local demux="$1"
    log "=== Running DADA2 denoising ==="
    
    # Set output variables
    OUT_DENOISE="$OUT_QIIME/rep-seqs-dada2.qza"
    OUT_TABLE="$OUT_QIIME/table-dada2.qza"
    OUT_TAXONOMY="$OUT_QIIME/taxonomy-dada2.qza"
    OUT_ALIGNED="$OUT_QIIME/aligned-rep-seqs-dada2.qza"
    OUT_MASKED="$OUT_QIIME/masked-aligned-rep-seqs-dada2.qza"
    OUT_TREE="$OUT_QIIME/fasttree-tree-dada2.qza"
    
    if [[ ! -f "$OUT_DENOISE" ]]; then
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$demux" \
            --p-trim-left-f 0 \
            --p-trunc-len-f 250 \
            --p-trim-left-r 0 \
            --p-trunc-len-r 250 \
            --p-n-threads "$THREADS" \
            --o-representative-sequences "$OUT_DENOISE" \
            --o-table "$OUT_TABLE" \
            --o-denoising-stats "$OUT_QIIME/dada2-stats.qza" \
            || error "DADA2 failed"
    fi
}

run_picrust2() {
    log "=== Running PICRUSt2 functional prediction with environment: $PICRUST_ENV ==="
    
    # Activate PICRUSt2 environment
    activate_env "$PICRUST_ENV"
    
    # Extract FASTA and BIOM
    local fasta_file=$(extract_from_qza "$OUT_DENOISE" "dna-sequences.fasta")
    local biom_file=$(extract_from_qza "$OUT_TABLE" "feature-table.biom")
    
    local picrust_out="$OUT_QIIME/picrust2_results"
    mkdir -p "$picrust_out"
    
    local original_dir=$(pwd)
    cd "$picrust_out" || error "Cannot cd to $picrust_out"
    
    # Run PICRUSt2
    if [[ ! -d "KO_metagenome_out" ]]; then
        picrust2_pipeline.py \
            -s "$fasta_file" \
            -i "$biom_file" \
            -o "$picrust_out" \
            -p "$THREADS" \
            || error "PICRUSt2 pipeline failed"
    fi
    
    # Add descriptions
    log "Adding descriptions to PICRUSt2 output..."
    if [[ -f "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" ]]; then
        add_descriptions.py \
            -i "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
            -m KO \
            -o "pathways_out/path_abun_unstrat_descripKO.tsv.gz" \
            || log "Warning: KO descriptions failed"
    fi
    
    if [[ -f "pathways_out/path_abun_unstrat.tsv.gz" ]]; then
        add_descriptions.py \
            -i "pathways_out/path_abun_unstrat.tsv.gz" \
            -m METACYC \
            -o "pathways_out/path_abun_unstrat_descrip.tsv.gz" \
            || log "Warning: METACYC descriptions failed"
    fi
    
    # Unzip description files
    for desc_file in "pathways_out/path_abun_unstrat_descripKO.tsv.gz" "pathways_out/path_abun_unstrat_descrip.tsv.gz"; do
        if [[ -f "$desc_file" ]]; then
            gunzip -k -f "$desc_file" || log "Warning: Failed to unzip $desc_file"
        fi
    done
    
    cd "$original_dir" || error "Cannot return to $original_dir"
    
    # Deactivate PICRUSt2 environment
    conda deactivate
    
    log "PICRUSt2 processing complete!"
}

#-----------------------------
# Main Script
#-----------------------------

# No arguments? Show help
if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

# Parse arguments
TEMP=$(getopt -o w:s:f:r:n:m:q:d:p:t:h \
    --long work-dir:,seqs-dir:,forward-pattern:,reverse-pattern:,n-char:,o-manifest:,o-qiime:,d-denoise:,p-picrust:,threads:,qiime-env:,picrust-env:,silva-dir:,silva-version:,forward-primer:,reverse-primer:,help \
    -- "$@")

if [[ $? -ne 0 ]]; then
    usage
    exit 1
fi

eval set -- "$TEMP"

while true; do
    case "$1" in
        -w|--work-dir)
            WDIR="$2"; shift 2 ;;
        -s|--seqs-dir)
            SDIR="$2"; shift 2 ;;
        -f|--forward-pattern)
            FORWARD_PATTERN="$2"; shift 2 ;;
        -r|--reverse-pattern)
            REVERSE_PATTERN="$2"; shift 2 ;;
        -n|--n-char)
            N_CHAR="$2"; shift 2 ;;
        -m|--o-manifest)
            OUT_MANI="$2"; shift 2 ;;
        -q|--o-qiime)
            OUT_QIIME="$2"; shift 2 ;;
        -d|--d-denoise)
            DENOISE="$2"; shift 2 ;;
        -p|--p-picrust)
            PICRUST="$2"; shift 2 ;;
        -t|--threads)
            THREADS="$2"; shift 2 ;;
        --qiime-env)
            QIIME_ENV="$2"; shift 2 ;;
        --picrust-env)
            PICRUST_ENV="$2"; shift 2 ;;
        --silva-dir)
            SILVA_DIR="$2"; shift 2 ;;
        --silva-version)
            SILVA_VERSION="$2"; shift 2 ;;
        --forward-primer)
            FORWARD_PRIMER="$2"; shift 2 ;;
        --reverse-primer)
            REVERSE_PRIMER="$2"; shift 2 ;;
        -h|--help)
            usage; exit 0 ;;
        --)
            shift; break ;;
        *)
            echo "Unknown option: $1" >&2
            usage; exit 1 ;;
    esac
done

# Validate arguments
validate_args

# Log start
log "=== Pipeline started ==="
log "Working directory: $WDIR"
log "Sequences directory: $SDIR"
log "Output directory: $OUT_QIIME"
log "Denoising method: $DENOISE"
log "PICRUSt2: $PICRUST"
log "Threads: $THREADS"
log "QIIME2 environment: $QIIME_ENV"
if [[ -n "$FORWARD_PRIMER" && -n "$REVERSE_PRIMER" ]]; then
    log "🔬 Primers provided: Forward='$FORWARD_PRIMER', Reverse='$REVERSE_PRIMER'"
    log "ℹ️  Will build region-specific classifier"
else
    log "⚠️  No primers provided"
    log "ℹ️  Will use full-length SILVA database (comprehensive but slower)"
fi
if [[ "$PICRUST" == "true" ]]; then
    log "PICRUSt2 environment: $PICRUST_ENV"
fi

# Check conda is available
check_command "conda"

# Run pipeline
create_manifest
setup_silva_database
run_qiime2

if [[ "$PICRUST" == "true" ]]; then
    run_picrust2#!/usr/bin/env bash
#===============================================================================
# Script: qiime2_picrust2_pipeline.sh
# Author: Tiago Mota
# Date: 2025-09-08
# Description:
#   QIIME2 16S rRNA processing pipeline with optional PICRUSt2 functional
#   prediction. Designed to run standalone or within SLURM jobs.
#
# Usage:
#   ./qiime2_picrust2_pipeline.sh [OPTIONS]
#
# Example SLURM submission script:
#   #!/bin/bash
#   #SBATCH --job-name=16S_pipeline
#   #SBATCH --cpus-per-task=10
#   #SBATCH --mem-per-cpu=10G
#   #SBATCH --time=168:00:00
#   
#   ./qiime2_picrust2_pipeline.sh \
#     --work-dir /project/microbiome \
#     --seqs-dir /project/data/reads \
#     --forward-pattern _R1 \
#     --reverse-pattern _R2 \
#     --n-char 5 \
#     --o-manifest results/manifest \
#     --o-qiime results/qiime2 \
#     --qiime-env qiime2-amplicon-2025.7 \
#     --picrust-env picrust2 \
#     --d-denoise dada2 \
#     --p-picrust true \
#     --threads 10
#===============================================================================

set -euo pipefail  # Exit on error, undefined variables, pipe failures

#-----------------------------
# Global Configuration
#-----------------------------
SCRIPT_NAME="$(basename "$0")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default values
DENOISE="deblur"
PICRUST="false"
THREADS=4
QIIME_ENV="qiime2-amplicon-2025.7"
PICRUST_ENV="picrust2"
SILVA_VERSION="138.2"
SILVA_TARGET="SSURef_NR99"
FORWARD_PRIMER="CCTACGGGNGGCWGCAG"
REVERSE_PRIMER="GGACTACNVGGGTWTCTAAT"
SILVA_DIR=""  # Will be set based on work-dir

#-----------------------------
# Functions
#-----------------------------
usage() {
    cat << EOF
Usage: ${SCRIPT_NAME} [OPTIONS]

Required:
  -w, --work-dir <path>          Working directory (absolute path)
  -s, --seqs-dir <path>          Directory containing FASTQ files
  -f, --forward-pattern <str>    Pattern for forward reads (e.g., _R1, _1.fq.gz)
  -r, --reverse-pattern <str>    Pattern for reverse reads (e.g., _R2, _2.fq.gz)
  -n, --n-char <int>             Characters to trim for sample ID (excludes patterns)
  -m, --o-manifest <path>        Output manifest file path (without extension)
  -q, --o-qiime <path>           Output directory for QIIME2 results

Optional:
  -d, --d-denoise <method>       Denoising method: deblur (default) or dada2
  -p, --p-picrust <bool>         Run PICRUSt2 functional prediction (true/false)
  -t, --threads <int>            Number of threads (default: 4)
  --qiime-env <name>             Conda environment name for QIIME2 (default: qiime2-amplicon-2025.7)
  --picrust-env <name>           Conda environment name for PICRUSt2 (default: picrust2)
  --silva-dir <path>             SILVA database directory (default: work-dir/../SILVA)
  --silva-version <ver>          SILVA version (default: 138.2)
  --forward-primer <seq>         Forward primer sequence (default: CCTACGGGNGGCWGCAG)
  --reverse-primer <seq>         Reverse primer sequence (default: GGACTACNVGGGTWTCTAAT)
  -h, --help                     Show this help message

Example:
  ${SCRIPT_NAME} \\
    --work-dir /project/microbiome \\
    --seqs-dir /project/data/reads \\
    --forward-pattern _R1 \\
    --reverse-pattern _R2 \\
    --n-char 5 \\
    --o-manifest results/manifest \\
    --o-qiime results/qiime2 \\
    --qiime-env qiime2-amplicon-2024.2 \\
    --picrust-env picrust2_v2.5 \\
    --d-denoise dada2 \\
    --p-picrust true \\
    --threads 8

SLURM Example:
  Create a separate SLURM submission script that calls this pipeline:
  
  #!/bin/bash
  #SBATCH --job-name=16S_pipeline
  #SBATCH --cpus-per-task=10
  #SBATCH --mem-per-cpu=10G
  #SBATCH --time=168:00:00
  
  ./${SCRIPT_NAME} --work-dir ... --threads \$SLURM_CPUS_PER_TASK

EOF
}

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

error() {
    echo "[ERROR] $*" >&2
    exit 1
}

check_command() {
    if ! command -v "$1" &> /dev/null; then
        error "Command '$1' not found. Please install or load the appropriate module."
    fi
}

check_conda_env() {
    local env_name="$1"
    if ! conda env list | grep -q "^${env_name} "; then
        error "Conda environment '$env_name' not found. Available environments:"
        conda env list
        return 1
    fi
    log "Conda environment '$env_name' found"
}

activate_env() {
    local env_name="$1"
    log "Activating conda environment: $env_name"
    
    # Initialize conda for bash
    if [[ -f "$CONDA_PREFIX/etc/profile.d/conda.sh" ]]; then
        source "$CONDA_PREFIX/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    else
        error "Cannot find conda.sh. Please ensure conda is installed and accessible."
    fi
    
    conda activate "$env_name" || error "Failed to activate environment: $env_name"
    log "Activated environment: $env_name (Python: $(python --version 2>&1))"
}

validate_args() {
    local ok=true

    # Check directories
    if [[ -z "$WDIR" ]]; then
        error "--work-dir is required"
    fi
    if [[ ! -d "$WDIR" ]]; then
        error "Directory '$WDIR' does not exist"
    fi

    if [[ -z "$SDIR" ]]; then
        error "--seqs-dir is required"
    fi
    if [[ ! -d "$SDIR" ]]; then
        error "Directory '$SDIR' does not exist"
    fi

    # Check patterns
    if [[ -z "$FORWARD_PATTERN" ]]; then
        error "--forward-pattern is required and cannot be empty"
    fi
    if [[ -z "$REVERSE_PATTERN" ]]; then
        error "--reverse-pattern is required and cannot be empty"
    fi

    # Check N_CHAR
    if [[ -z "$N_CHAR" ]]; then
        error "--n-char is required"
    fi
    if ! [[ "$N_CHAR" =~ ^[0-9]+$ ]]; then
        error "--n-char must be a non-negative integer (got '$N_CHAR')"
    fi

    # Check manifest output
    if [[ -z "$OUT_MANI" ]]; then
        error "--o-manifest is required and cannot be empty"
    else
        local outmdir=$(dirname "$OUT_MANI")
        mkdir -p "$outmdir" || error "Cannot create output directory '$outmdir'"
    fi

    # Check QIIME output
    if [[ -z "$OUT_QIIME" ]]; then
        error "--o-qiime is required and cannot be empty"
    else
        mkdir -p "$OUT_QIIME" || error "Cannot create QIIME output directory '$outqdir'"
    fi

    # Validate denoise method
    if [[ "$DENOISE" != "deblur" && "$DENOISE" != "dada2" ]]; then
        log "Warning: --d-denoise should be 'deblur' or 'dada2'. Using default: 'deblur'."
        DENOISE="deblur"
    fi

    # Validate PICRUST
    if [[ "$PICRUST" =~ ^(true|True|T|TRUE|false|False|F|FALSE)$ ]]; then
        if [[ "$PICRUST" =~ ^(true|True|T|TRUE)$ ]]; then
            PICRUST="true"
        else
            PICRUST="false"
        fi
    else
        log "Warning: --p-picrust should be true or false. Using default: 'false'."
        PICRUST="false"
    fi

    # Validate threads
    if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
        error "--threads must be a positive integer (got '$THREADS')"
    fi

    # Set SILVA_DIR if not provided
    if [[ -z "$SILVA_DIR" ]]; then
        SILVA_DIR="$WDIR/../SILVA"
    fi
    mkdir -p "$SILVA_DIR"

    return 0
}

setup_silva_database() {
    local classifier="${SILVA_DIR}/classifier.qza"
    local classifier_all="${SILVA_DIR}/classifier_all.qza"
    
    if [[ -f "$classifier" && -f "$classifier_all" ]]; then
        log "SILVA database already exists in $SILVA_DIR"
        return 0
    fi
    
    log "=== Downloading and building SILVA database (version ${SILVA_VERSION}) ==="
    local original_dir=$(pwd)
    cd "$SILVA_DIR" || error "Cannot cd to $SILVA_DIR"
    
    # Get Silva data
    if [[ ! -f "silva-${SILVA_VERSION}-ssu-nr99-seqs.qza" ]]; then
        qiime rescript get-silva-data \
            --p-version "${SILVA_VERSION}" \
            --p-target "${SILVA_TARGET}" \
            --p-include-species-labels \
            --o-silva-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs.qza" \
            --o-silva-taxonomy "silva-${SILVA_VERSION}-ssu-nr99-tax.qza" \
            || error "Failed to download SILVA data"
    fi
    
    # Clean sequences
    if [[ ! -f "silva-${SILVA_VERSION}-ssu-nr99-seqs-cleaned.qza" ]]; then
        qiime rescript cull-seqs \
            --i-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs.qza" \
            --o-clean-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs-cleaned.qza" \
            || error "Failed to clean sequences"
    fi
    
    # Filter unwanted taxa
    if [[ ! -f "silva-${SILVA_VERSION}-ssu-nr99-seqs-filt.qza" ]]; then
        qiime rescript filter-seqs-length-by-taxon \
            --i-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs-cleaned.qza" \
            --i-taxonomy "silva-${SILVA_VERSION}-ssu-nr99-tax.qza" \
            --p-labels Archaea Bacteria Eukaryota Mitochondria Chloroplast Unknown NA metagenome unidentified \
            --p-min-lens 900 1200 5000 5000 5000 5000 5000 5000 5000 \
            --o-filtered-seqs "silva-${SILVA_VERSION}-ssu-nr99-seqs-filt.qza" \
            --o-discarded-seqs "silva-${SILVA_VERSION}-ssu-nr99-seqs-discard.qza" \
            || error "Failed to filter sequences"
    fi
    
    # Dereplicate
    if [[ ! -f "silva-${SILVA_VERSION}-ssu-nr99-tax-derep-super.qza" ]]; then
        qiime rescript dereplicate \
            --i-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs-filt.qza" \
            --i-taxa "silva-${SILVA_VERSION}-ssu-nr99-tax.qza" \
            --p-mode 'uniq' \
            --p-threads "$THREADS" \
            --o-dereplicated-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs-derep-super.qza" \
            --o-dereplicated-taxa "silva-${SILVA_VERSION}-ssu-nr99-tax-derep-super.qza" \
            || error "Failed to dereplicate sequences"
    fi
    
    # Extract V3-V4 region
    if [[ ! -f "ref-seqs_silvaV3V4.qza" ]]; then
        qiime feature-classifier extract-reads \
            --i-sequences "silva-${SILVA_VERSION}-ssu-nr99-seqs-filt.qza" \
            --p-f-primer "$FORWARD_PRIMER" \
            --p-r-primer "$REVERSE_PRIMER" \
            --p-trunc-len 120 \
            --p-min-length 100 \
            --p-max-length 400 \
            --o-reads "ref-seqs_silvaV3V4.qza" \
            || error "Failed to extract V3-V4 region"
    fi
    
    # Build V3-V4 classifier
    if [[ ! -f "$classifier" ]]; then
        qiime feature-classifier fit-classifier-naive-bayes \
            --i-reference-reads "ref-seqs_silvaV3V4.qza" \
            --i-reference-taxonomy "silva-${SILVA_VERSION}-ssu-nr99-tax-derep-super.qza" \
            --o-classifier "$classifier" \
            || error "Failed to build V3-V4 classifier"
    fi
    
    # Build full-length classifier
    if [[ ! -f "$classifier_all" ]]; then
        qiime feature-classifier fit-classifier-naive-bayes \
            --i-reference-reads "silva-${SILVA_VERSION}-ssu-nr99-seqs-cleaned.qza" \
            --i-reference-taxonomy "silva-${SILVA_VERSION}-ssu-nr99-tax.qza" \
            --o-classifier "$classifier_all" \
            || error "Failed to build full-length classifier"
    fi
    
    cd "$original_dir" || error "Cannot return to $original_dir"
    log "SILVA database setup complete!"
}

create_manifest() {
    log "=== Creating manifest file ==="
    cd "$WDIR" || error "Cannot cd to $WDIR"
    
    # Get file lists
    local file_list=($SDIR/*)
    local for_seqs=($(printf '%s\n' "${file_list[@]}" | grep "$FORWARD_PATTERN" || true))
    local rev_seqs=($(printf '%s\n' "${file_list[@]}" | grep "$REVERSE_PATTERN" || true))
    
    if [[ ${#for_seqs[@]} -eq 0 ]] || [[ ${#rev_seqs[@]} -eq 0 ]]; then
        error "No files found matching patterns: forward='$FORWARD_PATTERN', reverse='$REVERSE_PATTERN'"
    fi
    
    if [[ ${#for_seqs[@]} -ne ${#rev_seqs[@]} ]]; then
        error "Number of forward (${#for_seqs[@]}) and reverse (${#rev_seqs[@]}) files do not match."
    fi
    
    # Extract replicate IDs
    local rep_id=()
    for f in "${for_seqs[@]}"; do
        local id=$(basename "$f" | sed -e "s/${FORWARD_PATTERN}//")
        rep_id+=("$id")
    done
    
    # Write manifest
    {
        echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"
        for ((i=0; i<${#rep_id[@]}; i++)); do
            echo -e "${rep_id[i]}\t${for_seqs[i]}\t${rev_seqs[i]}"
        done
    } > "${OUT_MANI}.tsv"
    
    log "Manifest created: ${OUT_MANI}.tsv (${#rep_id[@]} samples)"
}

extract_from_qza() {
    local qza_file="$1"
    local file_pattern="$2"
    local output_var_name="$3"
    
    if [[ ! -f "$qza_file" ]]; then
        error "QZA file not found: $qza_file"
    fi
    
    local folder_name=$(unzip -Z -1 "$qza_file" 2>/dev/null | head -1 | sed 's/\/[^/]*$//')
    if [[ -z "$folder_name" ]]; then
        error "Cannot determine folder structure in $qza_file"
    fi
    
    local extracted_file="$OUT_QIIME/${folder_name}/data/${file_pattern}"
    if [[ ! -f "$extracted_file" ]]; then
        log "Extracting $qza_file..."
        unzip -q "$qza_file" -d "$OUT_QIIME/" || error "Failed to extract $qza_file"
    fi
    
    echo "$extracted_file"
}

get_optimal_trunc_length() {
    local qzv_file="$1"
    local folder_name=$(unzip -Z -1 "$qzv_file" 2>/dev/null | head -1 | sed 's/\/[^/]*$//')
    
    if [[ -z "$folder_name" ]]; then
        error "Cannot determine folder structure in $qzv_file"
    fi
    
    local tsv_file="$OUT_QIIME/${folder_name}/data/forward-seven-number-summaries.tsv"
    if [[ ! -f "$tsv_file" ]]; then
        unzip -q "$qzv_file" -d "$OUT_QIIME/" || error "Failed to extract $qzv_file"
    fi
    
    local target_col=$(awk -F'\t' 'NR==2 {for(i=2;i<=NF;i++) if($i+0 < 9000) {print i-2; exit}}' "$tsv_file")
    target_col=${target_col:-0}
    
    echo "$target_col"
}

run_qiime2() {
    log "=== Starting QIIME2 processing with environment: $QIIME_ENV ==="
    
    # Activate QIIME2 environment
    activate_env "$QIIME_ENV"
    
    # Import sequences
    local demux="$OUT_QIIME/paired-end-demux.qza"
    if [[ ! -f "$demux" ]]; then
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path "${OUT_MANI}.tsv" \
            --output-path "$demux" \
            --input-format PairedEndFastqManifestPhred33V2 \
            || error "Failed to import sequences"
    fi
    
    # Denoising
    if [[ "$DENOISE" == "deblur" ]]; then
        run_deblur "$demux"
    else
        run_dada2 "$demux"
    fi
    
    # Taxonomy classification
    log "=== Classifying taxonomy ==="
    local classifier="$SILVA_DIR/classifier_all.qza"
    if [[ ! -f "$classifier" ]]; then
        error "Classifier not found: $classifier"
    fi
    
    qiime feature-classifier classify-sklearn \
        --i-classifier "$classifier" \
        --i-reads "$OUT_DENOISE" \
        --p-confidence 0.97 \
        --p-n-jobs "$THREADS" \
        --o-classification "$OUT_TAXONOMY" \
        || error "Taxonomy classification failed"
    
    # Phylogenetic tree
    log "=== Building phylogenetic tree ==="
    if [[ ! -f "$OUT_ALIGNED" ]]; then
        qiime alignment mafft \
            --i-sequences "$OUT_DENOISE" \
            --o-alignment "$OUT_ALIGNED" \
            --p-n-threads "$THREADS" \
            || error "MAFFT alignment failed"
    fi
    
    if [[ ! -f "$OUT_MASKED" ]]; then
        qiime alignment mask \
            --i-alignment "$OUT_ALIGNED" \
            --o-masked-alignment "$OUT_MASKED" \
            || error "Masking failed"
    fi
    
    if [[ ! -f "$OUT_TREE" ]]; then
        qiime phylogeny fasttree \
            --i-alignment "$OUT_MASKED" \
            --o-tree "$OUT_TREE" \
            --p-n-threads "$THREADS" \
            || error "FastTree failed"
    fi
    
    # Deactivate QIIME2 environment
    conda deactivate
    
    log "QIIME2 processing complete!"
}

run_deblur() {
    local demux="$1"
    log "=== Running Deblur denoising ==="
    
    # Merge paired ends
    local joined="$OUT_QIIME/demux-joined.qza"
    if [[ ! -f "$joined" ]]; then
        qiime vsearch merge-pairs \
            --i-demultiplexed-seqs "$demux" \
            --o-merged-sequences "$joined" \
            --o-unmerged-sequences "$OUT_QIIME/demux-unmerged.qza" \
            || error "VSEARCH merge failed"
    fi
    
    # Quality filter
    local filtered="$OUT_QIIME/demux-joined-filtered.qza"
    if [[ ! -f "$filtered" ]]; then
        qiime quality-filter q-score \
            --i-demux "$joined" \
            --p-min-quality 30 \
            --o-filtered-sequences "$filtered" \
            --o-filter-stats "$OUT_QIIME/demux-joined-filter-stats.qza" \
            || error "Quality filtering failed"
    fi
    
    # Visualize
    local qzv="$OUT_QIIME/demux-joined-filtered.qzv"
    if [[ ! -f "$qzv" ]]; then
        qiime demux summarize \
            --i-data "$filtered" \
            --o-visualization "$qzv" \
            || error "Demux summarization failed"
    fi
    
    # Get optimal truncation length
    local trim_length=$(get_optimal_trunc_length "$qzv")
    log "Optimal truncation length: $trim_length"
    
    # Set output variables
    OUT_DENOISE="$OUT_QIIME/rep-seqs-deblur.qza"
    OUT_TABLE="$OUT_QIIME/table-deblur.qza"
    OUT_TAXONOMY="$OUT_QIIME/taxonomy-silvaV3V4-deblur.qza"
    OUT_ALIGNED="$OUT_QIIME/aligned-rep-seqs-deblur.qza"
    OUT_MASKED="$OUT_QIIME/masked-aligned-rep-seqs-deblur.qza"
    OUT_TREE="$OUT_QIIME/fasttree-tree-deblur.qza"
    
    # Run deblur
    if [[ ! -f "$OUT_DENOISE" ]]; then
        qiime deblur denoise-16S \
            --i-demultiplexed-seqs "$filtered" \
            --p-trim-length "$trim_length" \
            --p-sample-stats \
            --p-jobs-to-start "$THREADS" \
            --o-representative-sequences "$OUT_DENOISE" \
            --o-table "$OUT_TABLE" \
            --o-stats "$OUT_QIIME/deblur-stats.qza" \
            || error "Deblur failed"
    fi
}

run_dada2() {
    local demux="$1"
    log "=== Running DADA2 denoising ==="
    
    # Set output variables
    OUT_DENOISE="$OUT_QIIME/rep-seqs-dada2.qza"
    OUT_TABLE="$OUT_QIIME/table-dada2.qza"
    OUT_TAXONOMY="$OUT_QIIME/taxonomy-silvaV3V4-dada2.qza"
    OUT_ALIGNED="$OUT_QIIME/aligned-rep-seqs-dada2.qza"
    OUT_MASKED="$OUT_QIIME/masked-aligned-rep-seqs-dada2.qza"
    OUT_TREE="$OUT_QIIME/fasttree-tree-dada2.qza"
    
    if [[ ! -f "$OUT_DENOISE" ]]; then
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$demux" \
            --p-trim-left-f 0 \
            --p-trunc-len-f 250 \
            --p-trim-left-r 0 \
            --p-trunc-len-r 250 \
            --p-n-threads "$THREADS" \
            --o-representative-sequences "$OUT_DENOISE" \
            --o-table "$OUT_TABLE" \
            --o-denoising-stats "$OUT_QIIME/dada2-stats.qza" \
            || error "DADA2 failed"
    fi
}

run_picrust2() {
    log "=== Running PICRUSt2 functional prediction with environment: $PICRUST_ENV ==="
    
    # Activate PICRUSt2 environment
    activate_env "$PICRUST_ENV"
    
    # Extract FASTA and BIOM
    local fasta_file=$(extract_from_qza "$OUT_DENOISE" "dna-sequences.fasta")
    local biom_file=$(extract_from_qza "$OUT_TABLE" "feature-table.biom")
    
    local picrust_out="$OUT_QIIME/picrust2_results"
    mkdir -p "$picrust_out"
    
    local original_dir=$(pwd)
    cd "$picrust_out" || error "Cannot cd to $picrust_out"
    
    # Run PICRUSt2
    if [[ ! -d "KO_metagenome_out" ]]; then
        picrust2_pipeline.py \
            -s "$fasta_file" \
            -i "$biom_file" \
            -o "$picrust_out" \
            -p "$THREADS" \
            || error "PICRUSt2 pipeline failed"
    fi
    
    # Add descriptions
    log "Adding descriptions to PICRUSt2 output..."
    if [[ -f "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" ]]; then
        add_descriptions.py \
            -i "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
            -m KO \
            -o "pathways_out/path_abun_unstrat_descripKO.tsv.gz" \
            || log "Warning: KO descriptions failed"
    fi
    
    if [[ -f "pathways_out/path_abun_unstrat.tsv.gz" ]]; then
        add_descriptions.py \
            -i "pathways_out/path_abun_unstrat.tsv.gz" \
            -m METACYC \
            -o "pathways_out/path_abun_unstrat_descrip.tsv.gz" \
            || log "Warning: METACYC descriptions failed"
    fi
    
    # Unzip description files
    for desc_file in "pathways_out/path_abun_unstrat_descripKO.tsv.gz" "pathways_out/path_abun_unstrat_descrip.tsv.gz"; do
        if [[ -f "$desc_file" ]]; then
            gunzip -k -f "$desc_file" || log "Warning: Failed to unzip $desc_file"
        fi
    done
    
    cd "$original_dir" || error "Cannot return to $original_dir"
    
    # Deactivate PICRUSt2 environment
    conda deactivate
    
    log "PICRUSt2 processing complete!"
}

#-----------------------------
# Main Script
#-----------------------------

# No arguments? Show help
if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

# Parse arguments
TEMP=$(getopt -o w:s:f:r:n:m:q:d:p:t:h \
    --long work-dir:,seqs-dir:,forward-pattern:,reverse-pattern:,n-char:,o-manifest:,o-qiime:,d-denoise:,p-picrust:,threads:,qiime-env:,picrust-env:,silva-dir:,silva-version:,forward-primer:,reverse-primer:,help \
    -- "$@")

if [[ $? -ne 0 ]]; then
    usage
    exit 1
fi

eval set -- "$TEMP"

while true; do
    case "$1" in
        -w|--work-dir)
            WDIR="$2"; shift 2 ;;
        -s|--seqs-dir)
            SDIR="$2"; shift 2 ;;
        -f|--forward-pattern)
            FORWARD_PATTERN="$2"; shift 2 ;;
        -r|--reverse-pattern)
            REVERSE_PATTERN="$2"; shift 2 ;;
        -n|--n-char)
            N_CHAR="$2"; shift 2 ;;
        -m|--o-manifest)
            OUT_MANI="$2"; shift 2 ;;
        -q|--o-qiime)
            OUT_QIIME="$2"; shift 2 ;;
        -d|--d-denoise)
            DENOISE="$2"; shift 2 ;;
        -p|--p-picrust)
            PICRUST="$2"; shift 2 ;;
        -t|--threads)
            THREADS="$2"; shift 2 ;;
        --qiime-env)
            QIIME_ENV="$2"; shift 2 ;;
        --picrust-env)
            PICRUST_ENV="$2"; shift 2 ;;
        --silva-dir)
            SILVA_DIR="$2"; shift 2 ;;
        --silva-version)
            SILVA_VERSION="$2"; shift 2 ;;
        --forward-primer)
            FORWARD_PRIMER="$2"; shift 2 ;;
        --reverse-primer)
            REVERSE_PRIMER="$2"; shift 2 ;;
        -h|--help)
            usage; exit 0 ;;
        --)
            shift; break ;;
        *)
            echo "Unknown option: $1" >&2
            usage; exit 1 ;;
    esac
done

# Validate arguments
validate_args

# Log start
log "=== Pipeline started ==="
log "Working directory: $WDIR"
log "Sequences directory: $SDIR"
log "Output directory: $OUT_QIIME"
log "Denoising method: $DENOISE"
log "PICRUSt2: $PICRUST"
log "Threads: $THREADS"
log "QIIME2 environment: $QIIME_ENV"
if [[ "$PICRUST" == "true" ]]; then
    log "PICRUSt2 environment: $PICRUST_ENV"
fi

# Check conda is available
check_command "conda"

# Run pipeline
create_manifest
setup_silva_database
run_qiime2

if [[ "$PICRUST" == "true" ]]; then
    run_picrust2
fi

# Clean up temporary files (optional)
log "Cleaning up temporary extracted files..."
find "$OUT_QIIME" -type d -name "data" -exec rm -rf {} \; 2>/dev/null || true

log "=== Pipeline completed successfully at $(date) ==="
log "Output directory: $OUT_QIIME"
