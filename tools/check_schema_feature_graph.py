#!/usr/bin/env python3
"""Keep native libraries outside portable builds and libxml outside SQLite-only builds."""
import os
import subprocess

xml_dependencies = {"libxml", "bindgen", "clang-sys"}
sqlite_dependencies = {"rusqlite", "libsqlite3-sys"}
for flags in [[], ["--no-default-features"], ["--no-default-features", "--features", "sqlite"]]:
    result = subprocess.run(
        ["cargo", "tree", "--locked", "--edges", "normal,build", "--prefix", "none", *flags],
        check=True, capture_output=True, text=True,
    )
    sqlite_enabled = "sqlite" in flags
    # vcpkg is also a build dependency of bundled libsqlite3-sys; its presence
    # in the SQLite-only graph does not imply libxml selection.
    forbidden = xml_dependencies | (set() if sqlite_enabled else sqlite_dependencies | {"vcpkg"})
    selected = {line.split()[0] for line in result.stdout.splitlines() if line.split()}
    assert not (selected & forbidden), selected & forbidden
    if sqlite_enabled:
        assert sqlite_dependencies <= selected, sqlite_dependencies - selected
    environment = os.environ.copy()
    # The exact binding's build script fails on this nonexistent override if it
    # runs. These intentional invalid paths must not affect portable features.
    environment.update(LIBXML2="/openms-deliberately-missing-libxml2", LIBCLANG_PATH="/openms-deliberately-missing-libclang")
    subprocess.run(["cargo", "check", "--locked", "-j2", *flags], env=environment, check=True)
    print("Optional native dependency boundary verified:", flags or ["default"], flush=True)
