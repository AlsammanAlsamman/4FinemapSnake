Prompt for Copilot / code generator
- DO NOT START CODING until the user explicitly types exactly: start code.
- Goal: produce one independent Snakemake rule (rules/<name>.smk) that uses bioconfigme to load YAML config inside the rule if desired, and calls an external helper script in scripts/ to perform the work. The helper script MUST NEVER read or parse any YAML files — it must receive all inputs, outputs and parameters exclusively via command-line arguments, environment variables, or explicit file paths passed by the rule.
- The rule’s main declared output must be a .done marker file (e.g., results/.../<name>.done) to indicate success.
- All variables required by the SNakemake rule (wildcards, config keys it expects to find via bioconfigme, file patterns, param names) must be declared inside the .smk file.
- Any scripts (bash, Python, R, etc.) must live in scripts/. The Snakemake rule may read YAML with bioconfigme and translate its values into CLI arguments to pass to the script, but the script must treat those arguments only as primitive types (strings, ints, paths) and must not itself open/parse YAML.
- The rule must log stdout/stderr to results/log/ and declare logs through the log directive or outputs.
- The rule must declare resources with default values: resources: mem_mb=32000, time="00:30:00", cores=2.
- Provide a scripts/submit_<rule_name>.sh that runs snakemake to execute only this rule (or target .done) so the rule can be tested independently.
- The submit script should run snakemake in local mode by default unless the user requests a cluster/profile mode (ask that as a yes/no).
- The rule and scripts must be self-contained so the helper script can be run directly from the command line with explicit args (i.e., show an example invocation that runs the helper script independently of Snakemake).
- The helper script must:
  - Validate required CLI args and exit non-zero with clear messages if arguments are missing or invalid.
  - If invoked with no args or with --help, print a usage summary of required arguments, their types, example values, and a full example command.
  - Never read YAML/Config files. If any structured data is required at runtime, it must be supplied by the rule as a file argument (e.g., a TSV/JSON/plain text) or as explicit CLI arguments — but again, the script itself must not parse YAML.
  - Minimize dependencies; prefer standard library. If third-party libraries are needed, list them and add them to requirements.txt (prompt must confirm whether to create/update requirements.txt).
  - Keep functions small and maintainable; separate parsing, validation, processing, and output-writing into small functions.
  - Produce the .done marker on success and write any output files to paths passed as arguments.
- Files to be created (Copilot should list these before writing any file contents; do NOT write file contents until the user says start code):
  - rules/<rule_name>.smk — the Snakemake rule, containing any bioconfigme usage and translating config keys into CLI args for the script. All wildcards and variables declared inside this file.
  - scripts/<helper_script>.(sh|py|R) — the helper program (language to be decided via yes/no question).
  - scripts/submit_<rule_name>.sh — runs snakemake to execute only this rule (.done target).
  - README-snippet or short instructions explaining how to run the script directly and how Snakemake will pass arguments.
  - Optionally requirements.txt if third-party libs are required for the helper script.
- Before coding behavior (must be enforced by the generator):
  - Ask the user at least five yes/no questions (the exact list below). Wait for the user’s answers. Show the assumptions and the exact CLI argument names and types the helper script will expect and the exact config keys the Snakemake rule will read via bioconfigme. Ask the user to confirm. Only after the user replies with exactly: start code should any file contents be written.
  - The generator must enumerate all assumptions about config keys and the minimal bioconfigme API it expects (for example: load_config(paths) -> dict, or bioconfigme.get('modules.tool')). The generator must ask the user to confirm or correct those assumptions before coding.
  - The generator must produce a before-coding checklist and ask the user to confirm.
- Required yes/no questions the generator must ask the user before coding (ask these exact questions and collect answers):
  1) Do you want the helper script written in Python by default? (yes/no)
  2) Should the helper script require all parameters to be passed as CLI arguments (i.e., fail if any required argument is missing)? (yes/no)
  3) Should the Snakemake rule be allowed to read/merge multiple YAML config files via bioconfigme and then pass primitive values to the helper script? (yes/no)
  4) Should the submit script run snakemake in local mode by default (yes) or use a cluster/profile wrapper (no = use cluster)? (yes/no)
  5) Should the helper script print a detailed usage and list the expected arguments when invoked with no args? (yes/no)
  6) Should the helper script accept environment variables in addition to CLI arguments for optional parameters? (yes/no)
  7) If third-party libraries are required for the helper script, should the generator create or update requirements.txt? (yes/no)
  8) Do you want strict wildcard validation and regex constraints in the Snakemake rule to avoid accidental wildcard matches? (yes/no)
  (At least five must be answered; include all above; the generator may ask additional questions if needed.)
- Before-coding checklist the generator must present to the user and confirm:
  - Exact helper script filename and language.
  - Exact Snakemake rule filename and the wildcard names it will declare.
  - Exact CLI arguments the helper script will accept: name, type, required/optional, default if any, and example values.
  - Exact config keys (YAML paths) the Snakemake rule will read using bioconfigme and how those keys map to CLI args for the script.
  - Paths for results/, logs/, scripts/, and submit script.
  - Whether requirements.txt will be created/updated.
  - Example independent invocation of the helper script (the generator will show this before coding).
- Assumptions that the generator must state before coding (and ask the user to confirm or correct):
  - Assumed bioconfigme API and return types (e.g., bioconfigme.load(paths) -> dict; bioconfigme.get(key, default=None)).
  - Assumed config YAML keys and names (e.g., modules.<tool>.image, analysis.samples, outputs.base_dir). List these keys explicitly and state which are required vs optional.
  - Assumption that Snakemake will pass file paths and parameter values as plain strings and that helper script will convert/cast types (int, float) as needed.
  - Assumption that the helper script should produce a .done marker file path provided as an argument.
  - Any default behaviors if optional config keys are missing (e.g., default memory, default cores).
- Output format the generator must use before writing code:
  1) Display the required yes/no questions above and prompt the user to answer them.
  2) Show a short list of the CLI arguments and YAML keys the generator assumes the rule and script will use (no code).
  3) Ask the user to respond with answers and then the single phrase "start code" when ready.
- Final strict note:
  - Under no circumstance should any file content or code be produced until the user answers the yes/no questions and types exactly: start code.
  - After the user provides answers and types start code, produce the described files (rules/<rule_name>.smk, scripts/<helper_script>, scripts/submit_<rule_name>.sh, README snippet, optional requirements.txt), ensuring the helper script never reads YAML and all parameters are taken from CLI args, env vars, or explicit files.
