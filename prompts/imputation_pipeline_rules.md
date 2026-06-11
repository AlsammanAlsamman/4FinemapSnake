# General Project Creation Rules and Guidelines

## Communication Rules
- **Do not write summaries**: Avoid explaining what you've done after completing tasks
- **README policy**: Only create main README.md once with brief description. Never update or create other README files unless explicitly requested
- **Wait for explicit permission**: Never start coding until user types exactly "start code"

## Project Structure Philosophy  
- **Simplicity first**: Keep structure simple and maintainable
- **Build incrementally**: Start with backbone, add features one by one
- **Independent components**: Each module/rule should work standalone
- **Clear organization**: Logical directory structure with clear purposes

## Configuration Management
- **bioconfigme centralization**: All configuration access must go through utils/bioconfigme.py
- **Two config files only**: analysis.yml and software.yml for maximum simplicity
- **No direct YAML reading**: Scripts receive all parameters via CLI arguments, never read YAML directly
- **Software module structure**: Each tool configured with module, version, and params blocks
- **Fail fast validation**: Missing or invalid config values should raise clear errors immediately

## Development Approach
- **Backbone first**: Create directory structure and basic configs before implementation
- **Step-by-step workflow**: Each rule should be independently executable with submit.sh
- **Sequential execution**: Clear pipeline steps with .done marker files for tracking
- **No documentation creation**: Do not create README files or documentation unless explicitly requested
- **User-driven features**: Only implement what user explicitly requests

## Code Organization Principles
- **Clear separation**: Configs, rules, scripts, utilities in separate directories
- **Independent rules**: Each Snakemake rule must be self-contained and runnable alone
- **bioconfigme integration**: All config access through centralized utility functions
- **CLI-only scripts**: Helper scripts accept parameters only via command-line arguments
- **Consistent naming**: Use descriptive, consistent file and directory names
- **Extensible design**: Structure should accommodate future additions easily

## Workflow Management
- **Independent execution**: Each rule can be executed individually via submit.sh
- **Sequential steps**: Clear pipeline progression: standardize → harmonize → extract → calculate → finemapping → report
- **Manual control first**: Step-by-step execution before any automation
- **Resource awareness**: Configure for SLURM environments with appropriate resource requests
- **Marker files**: Use .done files to track completion of each step
- **Error handling**: Build in appropriate error checking and logging

## User Interaction Rules
- **Ask clarifying questions**: Get requirements clear before starting
- **Present checklists**: Show what will be built before coding
- **Wait for confirmation**: Never assume user preferences
- **Respect constraints**: Honor user's specific requirements and limitations

## File Management
- **Preserve existing files**: Don't modify unless explicitly requested
- **Create new files thoughtfully**: Only create what's actually needed
- **Use appropriate tools**: Choose right tool for file creation/editing
- **Maintain consistency**: Follow established patterns within project

## Pipeline Execution Management
- **steps.sh file**: Central execution script containing sequential pipeline commands
- **Living documentation**: steps.sh serves as both execution guide and pipeline documentation
- **Update after each rule**: Add new commands to steps.sh when creating new rules
- **Test commands included**: steps.sh can contain commented testing commands for rule validation
- **Sequential execution**: Each step depends on previous step completion via .done marker files

## FineMapHub Specific Architecture
- **7-step pipeline**: Standardize → Soft Harmonization → Population Harmonization → Extract Loci → Calculate LD → Fine-mapping → Reporting
- **Multiple fine-mapping methods**: SuSiE-R, SuSiE-inf, FINEMAP, GCTA-COJO support
- **Population-specific analysis**: Support for EUR, AMR, EAS, AFR reference panels
- **Simplified config structure**: Replace complex configurations with analysis.yml and software.yml only
- **Target analysis format**: Each analysis defines gwas_table, population, snplist, region, and methods
- **steps.sh structure**: Contains commands like `./submit.sh --snakefile rules/rule_name.smk target.done`

## Results Directory Structure
- **Sequential organization**: Numbered directories (01_standardized, 02_harmonized, etc.) for clear pipeline progression
- **Method-specific subdirectories**: Fine-mapping methods separated (susieR/, susieinf/, finemap/, cojo/)
- **Organized logging**: Step-specific log directories for troubleshooting and monitoring
- **Standard structure**: 
  ```
  results/
  ├── 01_standardized/          # Standardized GWAS data with consistent column names
  ├── 02_harmonized/            # Harmonized data (soft harmonization with reference)
  ├── 03_popbased_harmonized/   # Population-based harmonization using reference panels
  ├── 04_loci/                  # Extracted loci regions around target SNPs
  ├── 05_ld_matrices/           # LD matrices for each locus and population
  ├── 06_finemapping/           # Fine-mapping results by method
  │   ├── susieR/, susieinf/, finemap/, cojo/
  ├── 07_reports/               # Comparison reports and visualization
  └── log/                      # Step-specific log directories
      ├── standardization/, harmonization/, popbased_harmonization/
      ├── loci_extraction/, ld_calculation/, finemapping/, reporting/
  ```