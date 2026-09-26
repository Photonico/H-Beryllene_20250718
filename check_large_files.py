#!/usr/bin/env python3
"""Find files larger than 99 MB and add unignored paths to .gitignore."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

THRESHOLD_BYTES = 99_000_000


def find_large_files(root: Path) -> list[Path]:
    large_files: list[Path] = []

    def raise_walk_error(error: OSError) -> None:
        raise error

    for directory, subdirectories, filenames in os.walk(
        root, topdown=True, onerror=raise_walk_error
    ):
        subdirectories[:] = [name for name in subdirectories if name != ".git"]
        current_directory = Path(directory)
        for filename in filenames:
            path = current_directory / filename
            if path.is_symlink():
                continue
            if path.stat().st_size > THRESHOLD_BYTES:
                large_files.append(path.relative_to(root))

    return large_files


def ignored_paths(root: Path, paths: list[Path]) -> set[str]:
    if not paths:
        return set()

    input_paths = b"\0".join(os.fsencode(path.as_posix()) for path in paths) + b"\0"
    result = subprocess.run(
        [
            "git",
            "-C",
            str(root),
            "check-ignore",
            "--no-index",
            "-z",
            "--stdin",
        ],
        input=input_paths,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode not in (0, 1):
        error_message = os.fsdecode(result.stderr).strip()
        raise RuntimeError(error_message or "git check-ignore failed")

    return {
        os.fsdecode(path)
        for path in result.stdout.split(b"\0")
        if path
    }


def add_to_gitignore(root: Path, paths: list[Path]) -> None:
    gitignore = root / ".gitignore"
    existing = gitignore.read_text(encoding="utf-8") if gitignore.exists() else ""
    with gitignore.open("a", encoding="utf-8", newline="") as file:
        if existing and not existing.endswith(("\n", "\r")):
            file.write("\n")
        for path in paths:
            file.write(f"/{path.as_posix()}\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Find files larger than 99 MB and add unignored paths to .gitignore."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="report missing ignore entries without modifying .gitignore",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    try:
        repository_root = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--show-toplevel"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        if Path(repository_root).resolve() != root:
            raise RuntimeError("The script must be located at the Git repository root.")

        large_files = find_large_files(root)
        ignored = ignored_paths(root, large_files)
        missing = [path for path in large_files if path.as_posix() not in ignored]

        print(f"Found {len(large_files)} file(s) larger than 99 MB.")
        for path in large_files:
            size_mb = (root / path).stat().st_size / 1_000_000
            status = "already ignored" if path.as_posix() in ignored else "not ignored"
            print(f"{size_mb:,.2f} MB\t{status}\t{path.as_posix()}")

        if missing:
            if args.dry_run:
                print(f"Dry run: {len(missing)} path(s) would be added to .gitignore.")
            else:
                add_to_gitignore(root, missing)
                print(f"Added {len(missing)} path(s) to .gitignore.")
        else:
            print("No new .gitignore entries needed.")
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
