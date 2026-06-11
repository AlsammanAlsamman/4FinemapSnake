
- DO NOT START CODING until the user explicitly types exactly: start code.
- Goal: produce a Snakemake project skeleton that uses utils/bioconfigme to load YAML config files from configs/. All config files are YAML (.yml). README.md is the only markdown file.

Project layout to create
- configs/
  - software.yml (YAML template for tool modules and parameters — exact content shown below)
  - analysis.yml (YAML template with results directory and analysis configuration — exact content shown below)
- scripts/
  - example helper script (template that accepts CLI args; must NOT read YAML files — it must get all values from arguments or env vars)
  - submit script (runs snakemake for the single rule)
- rules/
  - example_complete.smk (template rule that uses sys.path.append("utils") and bioconfigme functions)
- inputs/
  - (folder for SNP lists and input files - no README files)
- utils/
  - bioconfigme.py (module with documented API including get_results_dir() function)
- steps.sh (pipeline execution commands)
- submit.sh (main submission script)
- README.md (top-level, brief description only - never update unless explicitly requested)

Exact config templates (must be created verbatim)

1) configs/software.yml (exact YAML template)
# Software modules and tool configurations
plink2:
  module: "plink2/2.00a3.3lm"
  version: "2.00a3.3lm"
  params:
    memory: "--memory 32000"
    threads: "--threads 2"

r:
  module: "R/4.4.1-mkl"
  version: "4.4.1"
  params:
    susieR_max_causal: "10"
    coverage: "0.95"

python:
  module: "python/3.9"
  version: "3.9"

bcftools:
  module: "bcftools/1.15"
  version: "1.15"
  params:
    output_type: "-O z"
    threads: "--threads 2"

(Keep structure: tool_name -> module, version, params blocks)

2) configs/analysis.yml (exact YAML template)
# Analysis configuration and paths
project_name: "finemapping-analysis"
version: "0.1.0"
results_dir: "/path/to/results"  # Absolute path to results directory

# Reference files for harmonization
harmonization:
  ref_fasta: "/path/to/human_g1k_v37.fasta"
  ref_vcf: "/path/to/1000g.vcf.gz"

# Population reference panels
ref_panel:
  EUR: "/path/to/EUR_plink"
  AMR: "/path/to/AMR_plink"
  EAS: "/path/to/EAS_plink"

# GWAS datasets
gwastables:
  - name: example_gwas
    samples: 10000
    file: "/path/to/gwas.tsv"

# Analysis targets
target_analysis:
  example_analysis:
    gwas_table: example_gwas
    snplist: "inputs/snps.txt"
    population: EUR
    region: 1000000
    type: ["susieR", "susieinf"]

# Default resources
default_resources:
  mem_mb: 32000
  cores: 2
  time: "00:30:00"

(utils/bioconfigme must provide get_results_dir() function to access results_dir key)

Utilities and behavior requirements
- utils/bioconfigme must provide documented API including get_results_dir(), get_software_module(tool), get_analysis_value(path)
- All Snakemake rules must import bioconfigme using: sys.path.append("utils") then from bioconfigme import function_name
- Helper scripts must never open or parse any YAML file. They must accept all required inputs/params via CLI args, environment variables, or explicit data files passed as arguments.
- The example rule must:
  - be independent and live in rules/example_complete.smk
  - use sys.path.append("utils") at top before importing bioconfigme
  - declare explicit wildcards it expects (e.g., sample)
  - declare resources: mem_mb=32000, time="00:30:00", cores=2
  - declare logs under get_results_dir()/log/ path from bioconfigme
  - have primary output as a .done marker file
  - call an external script in scripts/ and pass all parameters as CLI args or env vars
- The submit script (scripts/submit_example.sh) must run snakemake targeting only the rule’s .done target; default to local mode unless user requests cluster/profile.
- If any code requires third-party libraries (e.g., PyYAML) the generator must add those to requirements.txt and list them explicitly.

Files to produce (generator must list before writing content; do NOT write files until user says "start code")
- configs/software.yml
- configs/analysis.yml
- utils/bioconfigme.py (with get_results_dir() function)
- rules/example_complete.smk (with sys.path.append("utils") import pattern)
- scripts/example_helper.py (or chosen language)
- scripts/submit_example.sh
- inputs/ (directory only - no README files)
- steps.sh (pipeline execution commands)
- submit.sh (main submission script)
- README.md (top-level, brief description - never update unless explicitly requested)

Assumptions the generator must state and ask the user to confirm before coding
- bioconfigme API: get_results_dir() -> str, get_software_module(tool) -> str, get_analysis_value(path) -> value
- bioconfigme import pattern: sys.path.append("utils") then from bioconfigme import function_name
- analysis.yml contains results_dir key with absolute path to results directory  
- software.yml structure: tool_name -> module, version, params blocks as shown
- Helper scripts must produce the .done marker at a path passed as an argument
- Default resources if missing in analysis.yml: mem_mb=32000, cores=2, time="00:30:00"
- inputs/ folder exists for SNP lists and input data files

Before-coding checklist the generator must present and confirm with the user
- Confirm helper script language (Python, Bash, R).
- Confirm helper script filename and exact CLI arguments (names/types/required).
- Confirm rule filename and wildcard names.
- Confirm paths for results/, logs/, scripts/, and configs/.
- Confirm whether requirements.txt should be created/updated for PyYAML (or other libs).
- Confirm submit script mode (local vs cluster/profile).

Output format required before writing code
1) Display the required yes/no candidate clarifying questions (see above) and prompt the user to answer them.
2) Show a short list of the CLI arguments and YAML keys the generator will expect (no code).
3) Ask the user to respond with answers and then the single phrase "start code" when ready.

Final notes
- All config files will be YAML (.yml) as requested; only README.md will be markdown.
- Helper scripts must never read YAML; Snakemake rules may read YAML via bioconfigme and must pass all primitive values to scripts as CLI args or env vars.
- I will not write any file contents until you answer the clarifying questions and type exactly: start code.