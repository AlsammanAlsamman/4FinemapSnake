# bioconfigme Architectural Baseline

## Purpose
- Centralize access to analysis and software configuration so workflows never hard-code paths, tool names, or parameters.
- Enforce consistent validation, logging, and error messaging when configuration values are missing or malformed.
- Provide a thin, reusable Python API that any Snakemake rule, helper script, or orchestration layer can import without pulling in unrelated project code.

## Design Principles
- **Single source of truth**: Load YAML (or other structured) configuration only through this module. No downstream component should read configuration files directly.
- **Fail fast**: Surface missing keys or invalid values immediately with actionable error messages. Prefer explicit `ValueError`/`RuntimeError` over silent defaults.
- **Composable helpers**: Expose focused functions that describe *intent* (e.g., `get_plink_conversion_params`) instead of returning raw config blobs.
- **Side-effect free**: Keep functions pure read-only utilities. Leave filesystem mutations, module loading, and command execution to callers.
- **Extensible**: New workflows should extend via additional helpers or by composing existing ones—never by forking the module.

## Public API Shape
- `get_results_dir() -> str`
  - Return the absolute path to results directory from analysis.yml results_dir key
  - Must be used for all output paths and log directories in rules
- `get_analysis_value(path: Sequence[str], *, required: bool = True, default: Optional[Any] = None) -> Any`
  - Resolve values from analysis configuration (e.g., `configs/analysis.yml`). Support nested lookups using path segments.
  - If `required=True`, raise with a context-rich message when the value is absent; otherwise return the provided default.
- `get_software_module(tool_name: str) -> str`
  - Return the environment/module declaration for a tool from software.yml
  - Encapsulate differences between cluster module systems, Conda, or container tags.
- Task-specific facades (e.g., `get_plink_conversion_params()`, `get_susie_params()`)
  - Wrap repetitive configuration lookups from software.yml tool params sections
  - Provide docstrings explaining their contract so downstream code knows which keys must exist.

## Usage Expectations
- Import pattern in Snakemake rules: `sys.path.append("utils")` then `from bioconfigme import get_results_dir, get_software_module`
- Always use get_results_dir() for output paths instead of hardcoding results directory
- Always destructure helper return values immediately for clarity, e.g.:
  ```python
  module_name, conversion_cfg = get_plink_conversion_params()
  ```
- When invoking from Bash, prefer tiny Python shims (or Snakemake `run` blocks) to avoid duplicating parsing logic in shell.
- Log the configuration values you consume (without leaking secrets) to aid reproducibility and debugging.

## Adding New Helpers
1. Identify repeated config access in the pipeline.
2. Define a dedicated helper that returns an immutable structure (tuple, `dataclass`, or `TypedDict`).
3. Validate inputs using `_require_config_value`-style utilities; raise meaningful errors.
4. Include minimal inline comments when logic is non-obvious (e.g., unit conversions or fallback precedence).
5. Update inline documentation in this file and relevant pipeline docs.

## Error Handling & Logging
- Use descriptive exception messages that mention the missing key path and the calling context.
- Do not print directly; leave logging to callers so they can route messages to Snakemake logs, SLURM output, or structured loggers.
- When returning compound objects, ensure downstream code can log them safely (e.g., convert `Path` to `str`).

## Testing & Validation
- Provide lightweight unit tests or doctests for new helpers whenever feasible.
- During development, add temporary assertions in consumers to confirm shapes and expected keys.
- For HPC environments, verify that `get_software_module` produces the exact string needed for `module load`/`ml` commands on the target cluster.

## Extension Playbook
- Need a new tool? Add software mapping to `configs/software.yml` with module, version, and params blocks, expose via `get_software_module`
- Need path configuration? Use `get_results_dir()` for results directory, extend with similar functions for other paths from analysis.yml
- Need scenario-specific defaults (e.g., dataset-specific overrides)? Extend the helper to merge base config while preserving backwards compatibility.
- Moving to a different config format? Isolate parsing behind a loader function so helper signatures remain stable.

## Anti-Patterns to Avoid
- Bypassing `bioconfigme` by directly reading YAML in scripts.
- Returning raw dictionaries without validation.
- Embedding filesystem or network side effects in helper functions.
- Duplicating helper logic across multiple modules.
- Hiding errors by returning `None` for required values.

## Documentation & Prompts
- Reference this file when drafting AI prompts so agents understand expected patterns before editing workflow code.
- Keep guidance project-agnostic: describe responsibilities and conventions rather than specific dataset details.
- Do not create or update README files unless explicitly requested by user.
- Only main README.md in root directory should exist with brief description - never update unless explicitly requested.
- When significant helpers or config schemas change, update this document first, then cascade to prompts only.
