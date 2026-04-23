#!/usr/bin/env python3
import argparse
import os


def main() -> None:
    parser = argparse.ArgumentParser(description="Example helper script that only consumes CLI args")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--module", required=True)
    parser.add_argument("--input-file", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.done_file), exist_ok=True)

    # No YAML reads here: everything is passed from the rule.
    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write(f"sample={args.sample}\n")
        handle.write(f"module={args.module}\n")
        handle.write(f"input={args.input_file}\n")


if __name__ == "__main__":
    main()
