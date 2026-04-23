# Pipeline Ancestry Guardrail Prompt

Use this prompt whenever modifying rules, scripts, or config related to ancestry processing.

You are updating a local-ancestry + Tractor pipeline that must support ANY number of ancestries.

Hard requirements:
1. Never hardcode ancestry count, labels, codes, color maps, or population names.
2. Always derive ancestry definitions from config at runtime.
3. Ensure all loops, matrix shapes, model inputs, and outputs scale with N ancestries.
4. Verify mapping consistency end-to-end:
   - population -> label
   - label -> numeric code (if used)
   - label -> color
   - code -> label (if used)
5. Before finalizing changes, run a mapping integrity check:
   - no missing keys
   - no duplicate labels
   - no duplicate numeric codes
   - every plotted/reported ancestry has a valid color
6. Ensure generated filenames and merged outputs do not assume fixed ancestries (no anc0..anc3 assumptions unless derived from N).
7. Ensure documentation and code comments match actual config behavior.
8. Preserve backward compatibility where possible.

Validation checklist (must pass before completion):
- [ ] Zero hardcoded ancestry set assumptions (e.g., AFR/AMR/EUR/EAS fixed set).
- [ ] Zero hardcoded fixed ancestry count assumptions (e.g., 4 ancestries).
- [ ] All ancestry-aware steps consume config-derived definitions.
- [ ] Plot legends and color assignment are deterministic and complete.
- [ ] Statistical models drop one configured reference ancestry only when required.
- [ ] Any default fallback is explicit and documented.

When you submit changes, include a short "Ancestry Safety Report" containing:
- Source of ancestry definitions in config
- Where N ancestries are iterated dynamically
- Mapping integrity checks added/confirmed
- Any residual risks or assumptions
