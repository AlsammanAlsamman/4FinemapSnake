

Important: DO NOT WRITE ANY CODE OR FILE CONTENTS until I type exactly: start code.

Prompt (give this to the generator)
- Goal: produce one independent Snakemake rule (rules/<name>.smk) that uses utils/bioconfigme to read config/ YAML and calls external helper script(s) in scripts/. The rule's primary output must be a .done marker file. Resources defaults: mem_mb=32000, time="00:30:00", cores=2. Logs must go to results/log/. Scripts must never read YAML — all arguments passed from the rule.

Key Architecture Rules
- Import bioconfigme using: sys.path.append("utils") then from bioconfigme import get_results_dir, get_software_module
- NO inline Python code in rule shell sections - use bioconfigme functions instead  
- Rule calls bash script with CLI arguments (no temp files or parameter files)
- Bash scripts contain actual tool commands (PLINK, bcftools, etc.)
- Snakemake handles SLURM parallelization - no sbatch calls inside rules
- Each rule should be independent and look for previous rule outputs as input
- Use get_results_dir() from bioconfigme to get results directory path from analysis.yml
- Update steps.sh with simple format: comment for step name + command only (no echo/notes)
- Command format in steps.sh: ./submit.sh --snakefile rules/RULE_NAME.smk TARGET_FILE
- steps.sh serves as living pipeline documentation and execution guide
- Include commented test commands in steps.sh for rule validation when appropriate

Pre-coding requirements (must follow)
1. Before any code, ask at least 5 yes/no questions. Make these questions adaptive: derive most of them from what you understood about my request and the config/layout I described. Include at least some of the fixed topics below (naming, wildcards, submit mode, helper language, validation), but tailor wording and add follow-ups where needed.
2. Present a very short (1–2 line) statement of what you understood from my request.
3. Show 2–3 proposed approaches or minor variants (one-line each) and label one recommended.
4. List the exact files you will create (filenames only).
5. List the config keys and wildcard names you assume the rule will use (brief bullet list).
- List assumptions (1–4 bullets) including expected bioconfigme API (e.g., get_results_dir()->str; get_software_module(tool)->str).
7. Present a before-coding checklist for me to confirm (language for helper script(s), results/log/scripts paths, rule name pattern, .done pattern, local vs cluster submit).
8. Wait for my answers to the yes/no questions and my confirmation of the checklist. Only after I reply with exactly the phrase start code should you generate any files or code.

Required flexible yes/no questions (generator must include at least these topics, but adapt wording and add error- or project-specific follow-ups)

Behavior notes for generator
- The generator must adapt further follow-up yes/no questions based on my answers and on what it infers from my earlier messages; do not use a rigid fixed list.
- Do not create any documentation, summaries, or README files unless I explicitly request them.
- Only create main README.md once with brief description - never update unless explicitly requested.
- Enumerate ALL config keys you will read from config/ (names and types) before coding and ask me to confirm them.
- If you make assumptions, state them explicitly and ask me to confirm or correct them.
- When I say start code, produce only the agreed files and code (no extra files) and ensure:
  - rules/<rule_name>.smk contains sys.path.append("utils") and bioconfigme imports at top
  - rules/<rule_name>.smk uses get_results_dir() for output paths and logs
  - scripts/ contains helper script(s) that accept CLI args (never read YAML)
  - scripts/submit_<rule_name>.sh runs snakemake to produce only the .done target
  - steps.sh updated with new command following format: # Step Name\n./submit.sh --snakefile rules/RULE_NAME.smk TARGET_FILE
  - Logs go to get_results_dir()/log/ directory
  - Resources defaults are included in the rule

Output format before coding (must be exactly)
1) The adaptive yes/no questions (show them).
2) One short sentence: "I understand: <your brief understanding>".
3) 2–3 one-line approach options (one recommended).
4) Files to create (filenames only).
5) Assumed config keys & wildcard names (brief bullets).
6) Before-coding checklist (items to confirm).
7) Prompt: "Answer the questions (yes/no) and then type exactly: start code"

Strict final note
- Under no circumstance produce file contents or code before I answer the questions and type exactly: start code.