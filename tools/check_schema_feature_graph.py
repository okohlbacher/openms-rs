#!/usr/bin/env python3
"""Prove default and no-feature builds do not select libxml build dependencies."""
import os
import subprocess

for flags in [[], ["--no-default-features"]]:
    result = subprocess.run(
        ["cargo", "tree", "--locked", "--edges", "normal,build", "--prefix", "none", *flags],
        check=True, capture_output=True, text=True,
    )
    forbidden = {"libxml", "bindgen", "clang-sys", "vcpkg"}
    selected = {line.split()[0] for line in result.stdout.splitlines() if line.split()}
    assert not (selected & forbidden), selected & forbidden
    environment = os.environ.copy()
    # The exact binding's build script fails on this nonexistent override if it
    # runs. These intentional invalid paths must not affect portable features.
    environment.update(LIBXML2="/openms-deliberately-missing-libxml2", LIBCLANG_PATH="/openms-deliberately-missing-libclang")
    subprocess.run(["cargo", "check", "--locked", "-j2", *flags], env=environment, check=True)
    print("No libxml native build work:", flags or ["default"], flush=True)
