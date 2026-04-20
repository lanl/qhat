"""Extract relevant equations from the Childs et al paper."""

import PyPDF2
import re
import os
import glob

# Find the PDF file
pdf_path = "/vast/home/bkkrueger/.claude/projects/-vast-home-bkkrueger-qc4ws-lisdi-git-qhat-worktrees-bkk-commutator-qre/1b372140-31cc-4899-b3c0-378f66f8c5be/tool-results/webfetch-1776705006385-jfcs70.pdf"

if not os.path.exists(pdf_path):
    print(f"ERROR: Could not find PDF file at {pdf_path}")
    exit(1)
print(f"Reading PDF: {pdf_path}\n")

# Open the PDF
with open(pdf_path, 'rb') as file:
    reader = PyPDF2.PdfReader(file)
    num_pages = len(reader.pages)
    print(f"Total pages: {num_pages}\n")

    # Search for pages containing "Equation" or "eq." and numbers around 144-146
    relevant_pages = []

    for page_num in range(num_pages):
        page = reader.pages[page_num]
        text = page.extract_text()

        # Look for equation references
        if re.search(r'\b(Equation|Eq\.|eq\.)?\s*\(?1[34][0-9]\)?', text, re.IGNORECASE):
            relevant_pages.append((page_num, text))

        # Also look for specific terms
        if any(term in text.lower() for term in ['second-order', 'nested commutator', 'c_1', 'c_2']):
            if page_num not in [p[0] for p in relevant_pages]:
                relevant_pages.append((page_num, text))

    print("="*70)
    print("PAGES POTENTIALLY CONTAINING EQUATIONS 144-145")
    print("="*70)

    for page_num, text in relevant_pages[:10]:  # Limit to first 10 matches
        print(f"\n{'='*70}")
        print(f"PAGE {page_num + 1}")
        print('='*70)

        # Extract lines containing equation numbers or relevant terms
        lines = text.split('\n')
        for i, line in enumerate(lines):
            # Check if line contains equation numbers or relevant keywords
            if (re.search(r'\b(Equation|Eq\.)?\s*\(?1[34][0-9]\)?', line, re.IGNORECASE) or
                any(term in line.lower() for term in ['c_1', 'c_2', 'nested commutator', 'second-order'])):

                # Print context (a few lines before and after)
                start = max(0, i - 2)
                end = min(len(lines), i + 3)

                print(f"\n--- Context around line {i} ---")
                for j in range(start, end):
                    marker = ">>> " if j == i else "    "
                    print(f"{marker}{lines[j]}")

        print("\n" + "-"*70)
        print("FULL PAGE TEXT:")
        print("-"*70)
        print(text[:2000])  # Print first 2000 chars of page
        if len(text) > 2000:
            print("\n... (truncated)")

print("\n" + "="*70)
print("SEARCHING FOR SPECIFIC FORMULAS")
print("="*70)

# Now do a more targeted search
for page_num in range(num_pages):
    page = reader.pages[page_num]
    text = page.extract_text()

    # Look for sum notation and commutators
    if 'i<j' in text or 'i<k' in text or '∑' in text:
        if '[H' in text or 'commutator' in text.lower():
            print(f"\n{'='*70}")
            print(f"PAGE {page_num + 1}: Contains sum notation with commutators")
            print('='*70)
            print(text)
            print("\n")
