#!/usr/bin/env python3
"""Scan Rmd for all rendering-breaking patterns."""
import re, sys

with open(sys.argv[1], 'r', encoding='utf-8') as f:
    content = f.read()

print("=== SCANNING FOR CORRUPTED LATEX ===")

# Find extit{ without backslash before it
# Simple approach: find all 'extit{' and check preceding char
issues = []
idx = 0
while True:
    idx = content.find('extit{', idx)
    if idx == -1:
        break
    # Check character before 'extit{'
    if idx > 0 and content[idx-1] != '\\':
        # Show context
        start = max(0, idx-3)
        end = min(len(content), idx+25)
        ctx = content[start:end]
        issues.append(f"  pos {idx}: ...{repr(ctx)}...")
    idx += 1

# Same for texttt{
idx = 0
while True:
    idx = content.find('texttt{', idx)
    if idx == -1:
        break
    if idx > 0 and content[idx-1] != '\\':
        start = max(0, idx-3)
        end = min(len(content), idx+25)
        issues.append(f"  pos {idx}: ...{repr(content[start:end])}...")
    idx += 1

if issues:
    print(f"FOUND {len(issues)} CORRUPTED COMMANDS:")
    for i in issues:
        print(i)
else:
    print("ZERO corrupted commands found")

# Count tabs
tabs = content.count('\t')
print(f"\nTab characters: {tabs}")

# Count proper commands
textit_ok = content.count('\\textit{')
texttt_ok = content.count('\\texttt{')
print(f"Proper \\textit: {textit_ok}")
print(f"Proper \\texttt: {texttt_ok}")
print(f"Total proper: {textit_ok + texttt_ok}")
